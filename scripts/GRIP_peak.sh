#!/bin/bash

#Find peaks for GRIP-seq

#GRIP-seq pipeline
#Author: Hongijang Liu
#Email: hongjiang.liu@ucsf.edu

# Exit the script if an error happens
set -e

if [[ $# -eq 0 ]]; then
  echo "Error: ROOTDIR must be provided as an input argument."
  exit 1
fi

REPEAT=$(cat ./repeat.txt)

ROOTDIR="$1"
ID=${2:-"current"}  # Default ID is current

echo "ROOTDIR: ${ROOTDIR} ..."
echo "ID: ${ID} ..."

echo "Making dir for peak..."
mkdir -p ${ROOTDIR}/peak/${ID}

echo "Preparing data form clipper..."
for NAME in ${REPEAT};
do
echo "Preparing data form clipper for ${NAME}..."
zcat ${ROOTDIR}/clipper/${ID}/${NAME}.peak.bed.gz | awk -F"\t" '{if ($6=="+") print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6; else print $1"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}' > ${ROOTDIR}/clipper/${ID}/${NAME}.summit.bed
done

echo "Calculating mapping depth..."
for NAME in ${REPEAT};
do
echo "Calculating mapping depth for ${NAME}..."
samtools depth -l 10 ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam > ${ROOTDIR}/STAR/${ID}/${NAME}.10.coverage
done

echo "Calling peaks..."
for NAME in ${REPEAT};
do
echo "Calling peaks for ${NAME}..."
Rscript ../peaks/peak.R \
${ROOTDIR}/clipper/${ID}/${NAME}.summit.bed \
${ROOTDIR}/STAR/${ID}/${NAME}.10.coverage \
1.4 \
${ROOTDIR}/peak/${ID}/${NAME}.GRIP.peak.bed
done


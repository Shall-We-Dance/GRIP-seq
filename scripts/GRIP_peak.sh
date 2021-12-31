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


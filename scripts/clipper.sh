#!/bin/bash

#Clipper

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
ID=${2:-"current"}  # Default ID is current


echo "ROOTDIR: ${ROOTDIR} ..."
echo "ID: ${ID} ..."

echo "Making dir for clipper..."
mkdir -p ${ROOTDIR}/clipper/${ID}

echo "Running clipper..."
for NAME in ${REPEAT}; 
do 
clipper -b ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam \
-s hg19 \
--processors \
-o ${ROOTDIR}/clipper/${ID}/${NAME}.peak.bed \
--FDR 0.01 --poisson-cutoff 1e-10 --minreads 5 --binomial 0.01
done

echo "Finished!"
exit 0

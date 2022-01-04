#!/bin/bash

#metaplotR

#GRIP-seq pipeline
#Author: Hongijang Liu
#Email: hongjiang.liu@ucsf.edu

# Exit the script if an error happens
set -e

if [[ $# -eq 0 ]]; then
  echo "Error: DIR must be provided as an input argument."
  exit 1
fi
REPEAT=$(cat ./repeat.txt)
ROOTDIR="$1"
TOOLDIR="$2"
ID=${3:-"current"}  # Default ID is current

echo "ROOTDIR: ${ROOTDIR} ..."
echo "ID: ${ID} ..."
echo "Making dir: ${ROOTDIR}/metaPlotR/${ID}..."
mkdir -p ${ROOTDIR}/metaPlotR/${ID}

for NAME in ${REPEAT}
do
echo "Drawing metaPlot for ${NAME}..."
sort -k1,1 -k2,2n ${ROOTDIR}/peak/${ID}/${NAME}.GRIP.peak.bed > ${ROOTDIR}/metaPlotR/${ID}/${NAME}.GRIP.peak.sorted.bed
intersectBed -a ${ROOTDIR}/metaPlotR/${ID}/${NAME}.GRIP.peak.sorted.bed -b ${TOOLDIR}/hg19_annot.sorted.bed -sorted -wo -s > ${ROOTDIR}/metaPlotR/${ID}/${NAME}.annot_m6a.sorted.bed
perl ${TOOLDIR}/rel_and_abs_dist_calc.pl --bed ${ROOTDIR}/metaPlotR/${ID}/${NAME}.annot_m6a.sorted.bed --regions ${TOOLDIR}/region_sizes.txt > ${ROOTDIR}/metaPlotR/${ID}/${NAME}.m6a.dist.measures.txt
#out
Rscript ../metaPlot/metaPlotR.r  ${ROOTDIR}/metaPlotR/${ID}/${NAME}.m6a.dist.measures.txt ${ROOTDIR}/metaPlotR/${ID}/${NAME}.pdf
done
echo "Finished!"

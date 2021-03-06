#!/bin/bash

#MEME

#GRIP-seq pipeline
#Author: Hongijang Liu
#Email: hongjiang.liu@ucsf.edu

if [[ $# -eq 0 ]]; then
  echo "Error: ROOTDIR & GENOMEDIR must be provided as an input argument."
  exit 1
fi

REPEAT=$(cat ./repeat.txt)
ROOTDIR="$1"
GENOMEDIR="$2"
ID=${3:-"current"}

mkdir -p ${ROOTDIR}/meme

for NAME in ${REPEAT}
do 
echo "Preparing for ${NAME}"
cat ${ROOTDIR}/peak/${ID}/${NAME}.GRIP.peak.bed | awk -F"\t" '{if ($6=="-") print $1"\t"$2-10"\t"$3+10}' > ${ROOTDIR}/meme/${ID}/${NAME}.broad-.bed
cat ${ROOTDIR}/peak/${ID}/${NAME}.GRIP.peak.bed | awk -F"\t" '{if ($6=="+") print $1"\t"$2-10"\t"$3+10}' > ${ROOTDIR}/meme/${ID}/${NAME}.broad+.bed
bedtools getfasta -fi ${GENOMEDIR}/hg19/GRCh37.p13.genome.fa \
-bed ${ROOTDIR}/meme/${ID}/${NAME}.broad-.bed \
-fo ${ROOTDIR}/meme/${ID}/${NAME}.broad-.fa
bedtools getfasta -fi ${GENOMEDIR}/hg19/GRCh37.p13.genome.fa \
-bed ${ROOTDIR}/meme/${ID}/${NAME}.broad+.bed \
-fo ${ROOTDIR}/meme/${ID}/${NAME}.broad+.fa
seqkit seq ${ROOTDIR}/meme/${ID}/${NAME}.broad-.fa -r -p -t dna> ${ROOTDIR}/meme/${ID}/${NAME}.broad-.r.fa
cat ${ROOTDIR}/meme/${ID}/${NAME}.broad-.r.fa ${ROOTDIR}/meme/${ID}/${NAME}.broad+.fa > ${ROOTDIR}/meme/${ID}/${NAME}.broad.fa
echo "Running meme for ${NAME}"
meme ${ROOTDIR}/meme/${ID}/${NAME}.broad.fa -dna -brief 1500000 -o  ${ROOTDIR}/meme/${ID}/${NAME} -time 14400 -mod zoops -nmotifs 3 -minw 5 -maxw 7 -objfun classic  -markov_order 0 
echo "meme finished for ${NAME}"
done

ROOTDIR=.
for NAME in WS013 WS014 WS015 WS016
do
#call_summits
zcat ${ROOTDIR}/${NAME}.peak.bed.gz |  awk -F"\t" '{if ($6=="+") print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6; else print $1"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}' > ${ROOTDIR}/${NAME}.summit.bed
done



Rscript ../find_peak.R \
/home/liuhongjiang/projects/CLIP/WS_211115/clipper/211214_2/WS015.summit.bed \
/home/liuhongjiang/projects/CLIP/WS_211115/STAR/211214/WS015.positive.10.coverage \
1.4 \
/home/liuhongjiang/projects/CLIP/WS_211115/find_peak/WS015.1.4.summits.filtered.bed

ROOTDIR=.
for NAME in WS013 WS014 WS015 WS016
do
#call_summits
zcat ${ROOTDIR}/${NAME}.peak.bed.gz |  awk -F"\t" '{if ($6=="+") print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6; else print $1"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}' > ${ROOTDIR}/${NAME}.summit.bed
done

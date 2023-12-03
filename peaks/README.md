# GRIP-seq
## GRIP-seq Peaks Calling

+ **Call peaks using [Clipper](https://github.com/YeoLab/clipper).**

```
# activate conda evironment - clipper3
conda activate clipper3

clipper -b ${ROOTDIR}/STAR/${INID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam \
-s hg19 \
-o ${ROOTDIR}/clipper/${ID}/${NAME}.peak.bed \
--FDR 0.01 --poisson-cutoff 1e-50 --minreads 5 --binomial 0.01
```

+ **Separate summits from clipper results.**

```
zcat ${ROOTDIR}/clipper/${ID}/${NAME}.peak.bed.gz |  awk -F"\t" '{if ($6=="+") print $1"\t"$2+2"\t"$2+3"\t"$6; else print $1"\t"$3-3"\t"$3-2"\t"$6}' > ${ROOTDIR}/clipper/${ID}/${NAME}.summit.bed

```

+ **Call peaks for GRIP-seq using our R scirpts [peaks.R](peaks.R).**

```
#arg1 <- Input bed file separate summits form clipper results
#arg2 <- reference depth from STAR result
#arg3 <- threshhold valuse, 2 will works for most case
#arg4 <- Out put bed file, the peaks for GRIP-seq
Rscripts peaks.R ${ROOTDIR}/clipper/${ID}/${NAME}.summit.bed
```

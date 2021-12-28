# GRIP-seq
## call peaks for GRIP-seq

+ **Call peaks using [Clipper](https://github.com/YeoLab/clipper).**

```
# activate conda evironment - clipper3
conda activate clipper3

clipper -b ${ROOTDIR}/STAR/${INID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam \
-s hg19 \
-o ${ROOTDIR}/clipper/${ID}/${NAME}.peak.bed \
--FDR 0.01 --poisson-cutoff 1e-50 --minreads 5 --binomial 0.01
```

+ **Separate summits form clipper results.**

+ **Call peaks for GRIP-seq using our R scirpts [peaks.R](peaks.R).**

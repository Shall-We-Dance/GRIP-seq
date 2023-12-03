# GRIP-seq scripts

## generate_genome_index.sh

### Generate genome indexes for [STAR](https://github.com/alexdobin/STAR)

  ```
  #basic usage
  bash scripts/generate_genome_index.sh ${GENOME_DIR}
  ```
  This will generate a hg19 genome index using defalut settings, which uses `--sjdbOverhang=100` and runs on 8 threads.
  
  To specify the threads used to generate, run:
  
  ```
  #specify CPU threads
  bash scripts/generate_genome_index.sh ${GENOME_DIR} ${CPU_THREADS}
  ```
  
  The ideal value of `--sjdbOverhang` is `max(ReadLength)-1`, to specify, run:
  
  ```
  #specify CPU threads & index length
  bash scripts/generate_genome_index.sh ${GENOME_DIR} ${CPU_THREADS} ${INDEX_LENGTH}
  ```
  

### 3.  Make directory for GRIP-seq pipeline.
  
```
mkdir -p ${ANALYSIS_DIR}/fastp/${ID}
mkdir -p ${ANALYSIS_DIR}/STAR/${ID}
mkdir -p ${ANALYSIS_DIR}/bigWig/${ID}
mkdir -p ${ANALYSIS_DIR}/clipper/${ID}
mkdir -p ${ANALYSIS_DIR}/output/${ID}
```
  
  The directory structure should be like: 
  
```
${TOOLSDIR}/
    meme/
        #Using meme to find motif
    STAR/
        #Using STAR to map reads
    clipper/
        #Using clipper to call peaks
    metaPlotR/
        #Using clipper to create metagene plots
${ANALYSIS_DIR}/
    raw_data/
        # Your reads files
    fastp/
        ${ID}/
            #fastp is used to filter raw data
    STAR/
        ${ID}/
            #STAR is used to map reads
    bigWig/
        ${ID}/
            #convert the mapping result to bigwig format for viewing
    clipper/
        ${ID}/
            #clipper is used to call peaks
    metaPlotR/
        ${ID}/
    output/
        ${ID}/
            #peak reads, meme results
${GENOME_DIR}/
    STAR_index/
        #the genome index for STAR
    hg19/
        #hg19 will be used to map reads by STAR

```

### 4.  GRIP-seq pipeline

```
# activate conda evironment - GRIP-seq
conda activate GRIP-seq
```

+ **Remove PCR deduplication and cut adapter using [fastp](https://github.com/OpenGene/fastp).**

```

## Read length filter

fastp -i ${ROOTDIR}/raw_data/${NAME}_L003_R1_001.fastq.gz \
-I ${ROOTDIR}/raw_data/${NAME}_L003_R2_001.fastq.gz \
-o ${ROOTDIR}/fastp/${ID}/${NAME}_R1_length.fastq.gz \
-O ${ROOTDIR}/fastp/${ID}/${NAME}_R2_length.fastq.gz \
-h ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_length.html \
-j ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_length.json \
--thread ${THREAD} -Q -A \
--length_required 101 \
--length_limit 101

## quality control and cut illumina adapters

fastp -i ${ROOTDIR}/fastp/${ID}/${NAME}_R1_length.fastq.gz \
-I ${ROOTDIR}/fastp/${ID}/${NAME}_R2_length.fastq.gz \
-o ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_dedup_adapter.fastq.gz \
-O ${ROOTDIR}/fastp/${ID}/${NAME}_R2_fastp_dedup_adapter.fastq.gz \
-h ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_dedup_adapter.html \
-j ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_dedup_adapter.json \
--thread ${THREAD} \
--dedup

## cut GRIP-seq own adapters

fastp -i ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_dedup_adapter.fastq.gz \
-I ${ROOTDIR}/fastp/${ID}/${NAME}_R2_fastp_dedup_adapter.fastq.gz \
-o ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_final.fastq.gz \
-O ${ROOTDIR}/fastp/${ID}/${NAME}_R2_astp_final.fastq.gz \
-h ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_final.html \
-j ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_final.json \
--thread ${THREAD} \
-A \
-f 12 -t 10 -F 10 -T 12
```

+ **Map reads using [STAR](https://github.com/alexdobin/STAR).**

```
STAR  --runThreadN ${CPU_THREADS} \
--genomeDir  ${GENOMEDIR_STAR} \
--readFilesCommand zcat \
--readFilesIn  ${ROOTDIR}/cutadapt/${ID}/${NAME}_R2_cutadapt_fastp_dedup_adapter.fastq.gz \
--outFileNamePrefix ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bam \
--outSAMtype BAM SortedByCoordinate ;
```

  If you meet the `File size limit exceeded(core dumped)` error, run command below to fix it. (This usually happens if `${CPU_THREADS} > 16`)
  ```
  ulimit -n 65535
  ```

+ **Index the mapping output `.bam` file using [samtools](https://www.htslib.org).**

```
samtools index -b ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam;
```

+ **Generate `.bw` file (bigwig) for mapping output using [deepTools](https://github.com/deeptools/deepTools).**

```
bamCoverage -b ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam \
-o ${ROOTDIR}/bigWig/${ID}/${NAME}.bigwig;
```

+ Call peaks using [our protocol](/peaks).

  To find out the summits of the peaks, run:
  ```
  zcat ${ROOTDIR}/${NAME}.peak.bed.gz | awk -F"\t" '{if ($6=="+") print $1"\t"$2+2"\t"$2+3"\t"$4"\t"$5"\t"$6; else print $1"\t"$3-3"\t"$3-2"\t"$4"\t"$5"\t"$6}' > ${ROOTDIR}/${NAME}.summit.bed
  sort -k1,1 -k2,2n ${ROOTDIR}/${NAME}.summit.bed > ${ROOTDIR}/${NAME}.summit.sorted.bed
  ```
  
  📒NOTE: We do not recommend you to use [MACS](https://github.com/macs3-project/MACS) for peak calling, as our results significantly differ from MACS's model. I attempted to use MACS, but it yielded weak peaks.

# GRIP-seq

Here is the repositories of our unpublished paper, which introduciong a new technique **GRIP-seq**.

You can use the pipeline described below to analysis GRIP-seq data.

## Install

1. Clone the package
```
git clone https://github.com/Shall-We-Dance/GRIP-seq.git
cd GRIP-seq
```

2. Create conda environment using `GRIP-seq.yml` and `clipper3.yml`.
```
# create conda environment for GRIP-seq
conda env create -f GRIP-seq.yml

# create conda environment for clipper3
conda env create -f clipper3.yml

# activate conda evironment - GRIP-seq
conda activate GRIP-seq
```
## Pipeline

### Usage

1.  Make directories.

+ Make a directory `ANALYSIS_DIR` for this analysis, put your raw data at `ANALYSIS_DIR/raw_data` folder. 

+ Make a directory `GENOME_DIR` for [STAR](https://github.com/alexdobin/STAR), which is a important tool used in our analysis. 

  The structure should be like:
  
```
${TOOLSDIR}/
    meme/
        #Using meme to find motif
    STAR/
        #Using STAR to map reads
    clipper/
        #Using clipper to call peaks
${ANALYSIS_DIR}/
    raw_data/
        # Your reads files
${GENOME_DIR}/
    hg19/
        #hg19 will be used to map reads by STAR

```

2.  Generate genome indexes for [STAR](https://github.com/alexdobin/STAR)


3.  Make directory for GRIP-seq pipeline.
  
```
mkdir -p ${ANALYSIS_DIR}/fastp/${ANALYSIS_DIR}
mkdir -p ${ANALYSIS_DIR}/cutadapt/${ANALYSIS_DIR}
mkdir -p ${ANALYSIS_DIR}/STAR/${ANALYSIS_DIR}
mkdir -p ${ANALYSIS_DIR}/bigWig/${ANALYSIS_DIR}
mkdir -p ${ANALYSIS_DIR}/macs2/${ANALYSIS_DIR}
mkdir -p ${ANALYSIS_DIR}/output/${ANALYSIS_DIR}
```
  
2.  Remove PCR deduplication and cut adapter using [fastp](https://github.com/OpenGene/fastp).

```
fastp -i ${ROOTDIR}/raw_data/${ID}/${NAME}_R1.fastq.gz \
-I ${ROOTDIR}/raw_data/${ID}/${NAME}_R2.fastq.gz \
-o ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_dedup_adapter.fastq.gz \
-O ${ROOTDIR}/fastp/${ID}/${NAME}_R2_fastp_dedup_adapter.fastq.gz \
-h ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_dedup_adapter.html \
-j ${ROOTDIR}/fastp/${ID}/${NAME}_fastp_dedup_adapter.json \
--thread ${THREAD} \
--dedup
```

3.  Cut adapter using [cutadapt](https://github.com/marcelm/cutadapt).

```
cutadapt -j ${THREAD} \
-u 12 -u -10 \
-o ${ROOTDIR}/cutadapt/${ID}/${NAME}_R1_cutadapt_fastp_dedup_adapter.fastq.gz ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_dedup_adapter.fastq.gz;
cutadapt -j ${THREAD} \
u 10 -u -12 \
-o ${ROOTDIR}/cutadapt/${ID}/${NAME}_R2_cutadapt_fastp_dedup_adapter.fastq.gz ${ROOTDIR}/fastp/${ID}/${NAME}_R2_fastp_dedup_adapter.fastq.gz;
```

4.  Map reads using [STAR](https://github.com/alexdobin/STAR).

```
STAR  --runThreadN 40 \
--genomeDir  ${GENOMEDIR_STAR} \
--readFilesCommand zcat \
--readFilesIn  ${ROOTDIR}/cutadapt/${ID}/${NAME}_R2_cutadapt_fastp_dedup_adapter.fastq.gz \
--outFileNamePrefix ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bam \
--outSAMtype BAM SortedByCoordinate ;
```

5.  Index the mapping output `.bam` file using [samtools](https://www.htslib.org).

```
samtools index -b ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam;
```

6.  Generate `.bw` file (bigwig) for mapping output using [deepTools](https://github.com/deeptools/deepTools).

```
bamCoverage -b ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam \
-o ${ROOTDIR}/bigWig/${ID}/${NAME}.bigwig;
```

7.  Call peaks using [Clipper](https://github.com/YeoLab/clipper).

```

```
  
  📒NOTE: We do not recommend using [MACS](https://github.com/macs3-project/MACS) to call peaks, our result is quite different form the model used by MACS, I've tried that tools but the result is week.

8.  Intersect output `.bam` with RefSeq to identify mRNA reads using [bedtools](https://github.com/arq5x/bedtools2).

```

```


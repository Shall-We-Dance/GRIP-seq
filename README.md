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

1.  Make directories.

+ Make a directory `ANALYSIS_DIR` for this analysis, put your raw data at `ANALYSIS_DIR/raw_data` folder. 

+ Make a directory `GENOME_DIR` for [STAR](https://github.com/alexdobin/STAR) to generate genome index. 

+ Make a directory `TOOLS_DIR`, and install these 4 tools: [STAR](https://github.com/alexdobin/STAR), [meme](https://meme-suite.org/meme/doc/download.html), [clipper](https://github.com/YeoLab/clipper), [metaPlotR](https://github.com/olarerin/metaPlotR) . 

  The directory structure should be like:
  
```
${TOOLS_DIR}/
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
${GENOME_DIR}/

```

2.  Generate genome indexes for [STAR](https://github.com/alexdobin/STAR)

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

3.  Make directory for GRIP-seq pipeline.
  
```
mkdir -p ${ANALYSIS_DIR}/fastp/${ID}
mkdir -p ${ANALYSIS_DIR}/cutadapt/${ID}
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
    cutadapt/
        ${ID}/
            #cutadapt is used to cut random adapter
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

4.  GRIP-seq pipeline


```
# activate conda evironment - GRIP-seq
conda activate GRIP-seq
```

+ Remove PCR deduplication and cut adapter using [fastp](https://github.com/OpenGene/fastp).

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

+ Cut adapter using [cutadapt](https://github.com/marcelm/cutadapt).

```
cutadapt -j ${THREAD} \
-u 12 -u -10 \
-o ${ROOTDIR}/cutadapt/${ID}/${NAME}_R1_cutadapt_fastp_dedup_adapter.fastq.gz ${ROOTDIR}/fastp/${ID}/${NAME}_R1_fastp_dedup_adapter.fastq.gz;
cutadapt -j ${THREAD} \
u 10 -u -12 \
-o ${ROOTDIR}/cutadapt/${ID}/${NAME}_R2_cutadapt_fastp_dedup_adapter.fastq.gz ${ROOTDIR}/fastp/${ID}/${NAME}_R2_fastp_dedup_adapter.fastq.gz;
```

+ Map reads using [STAR](https://github.com/alexdobin/STAR).

```
STAR  --runThreadN ${CPU_THREADS} \
--genomeDir  ${GENOMEDIR_STAR} \
--readFilesCommand zcat \
--readFilesIn  ${ROOTDIR}/cutadapt/${ID}/${NAME}_R2_cutadapt_fastp_dedup_adapter.fastq.gz \
--outFileNamePrefix ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bam \
--outSAMtype BAM SortedByCoordinate ;
```

+ Index the mapping output `.bam` file using [samtools](https://www.htslib.org).

```
samtools index -b ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam;
```

+ Generate `.bw` file (bigwig) for mapping output using [deepTools](https://github.com/deeptools/deepTools).

```
bamCoverage -b ${ROOTDIR}/STAR/${ID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam \
-o ${ROOTDIR}/bigWig/${ID}/${NAME}.bigwig;
```

+ Call peaks using [Clipper](https://github.com/YeoLab/clipper).

```
# activate conda evironment - clipper3
conda activate clipper3

clipper -b ${ROOTDIR}/STAR/${INID}/${NAME}/${NAME}.bamAligned.sortedByCoord.out.bam \
-s hg19 \
-o ${ROOTDIR}/clipper/${ID}/${NAME}.peak.bed \
--FDR 0.01 --poisson-cutoff 1e-50 --minreads 5 --binomial 0.01
```
  
  📒NOTE: We do not recommend using [MACS](https://github.com/macs3-project/MACS) to call peaks, our result is quite different form the model used by MACS, I've tried that tools but the result is week.

## Making UpSetPlot

Using [UpSetR](https://github.com/hms-dbmi/UpSetR).

To finish this step, please refer to `UpSetR` [folder](https://github.com/Shall-We-Dance/GRIP-seq/tree/main/UpSetR).

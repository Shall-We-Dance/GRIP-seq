# GRIP-seq

Here is the repositories of our unpublished paper, which introduciong a new technique - **GRIP-seq**.

You can use the pipeline described below to analysis GRIP-seq data.

## Install

1. Clone the package
```
git clone https://github.com/Shall-We-Dance/GRIP-seq.git
cd GRIP-seq
```

2. Create conda environment using `GRIP-seq.yml` file.
```
# create conda environment for GRIP-seq
conda env create -f GRIP-seq.yml

# activate conda evironment - GRIP-seq
conda activate GRIP-seq
```
## Pipeline

### Usage

1.  Make a directory `<ANALYSIS_DIR>` for this analysis, put your raw data at `<ANALYSIS_DIR>/raw_data` folder. 

2.  Make a directory `<GENOME_DIR>` for [STAR](https://github.com/alexdobin/STAR), which is a important tool used in our analysis. 

  The structure should be like:

3.  To use our GRIP-seq pipeline, run `GRIP-seq.sh` with presets below.

<ANALYSIS_DIR>
<GENOME_DIR>
<READS_NAME>
<LENGTH>

```
bash GRIP-seq.sh --dir <ANALYSIS_DIR> --
```

```
###################
##GRIP-seq manual##
###################

Version: v1.0.0
Code: https://github.com/Shall-We-Dance/GRIP-seq

Usage:  bash GRIP-seq.sh [options]

The options include:
  --help  Print this help menu.
  --dir The root dir <ANALYSIS_DIR> of our analysis
```

### Steps

1.  Make directory for GRIP-seq pipeline.
  
```
mkdir -p <ANALYSIS_DIR>/fastp/<ANALYSISI_ID>
mkdir -p <ANALYSIS_DIR>/cutadapt/<ANALYSISI_ID>
mkdir -p <ANALYSIS_DIR>/STAR/<ANALYSISI_ID>
mkdir -p <ANALYSIS_DIR>/bigWig/<ANALYSISI_ID>
mkdir -p <ANALYSIS_DIR>/macs2/<ANALYSISI_ID>
mkdir -p <ANALYSIS_DIR>/output/<ANALYSISI_ID>
```
  
2.  Remove PCR deduplication using [fastp](https://github.com/OpenGene/fastp).

```
fastp -i <ANALYSIS_DIR>/raw_data/<READS_NAME>.fastq.gz -o <ANALYSIS_DIR>/fastp/<READS_NAME>_fastp_dedup.fastq.gz --dedup --length_required <LENGTH> --disable_adapter_trimming -Q -h <ANALYSIS_DIR>/fastp/<READS_NAME>_fastp_dedup.html -j <ANALYSIS_DIR>/fastp/<READS_NAME>_fastp_dedup.json
```

3.  Cut adapter using [cutadapt](https://github.com/marcelm/cutadapt).

```

```

4.  Map reads using [STAR](https://github.com/alexdobin/STAR).

```

```

5.  Index the mapping output `.bam` file using [samtools](https://www.htslib.org).

```

```

6.  Generate `.bw` file (bigwig) for mapping output using [deepTools](https://github.com/deeptools/deepTools).

```

```

7.  Call peaks using [MACS](https://github.com/macs3-project/MACS).

```

```

8.  Intersect output `.bam` with RefSeq to identify mRNA reads using [bedtools](https://github.com/arq5x/bedtools2).

```

```


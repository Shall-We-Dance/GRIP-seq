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

To use our GRIP-seq pipeline, run `GRIP-seq.sh` with presets below.

```
###################
##GRIP-seq manual##
###################
Version: v1.0.0
Code: https://github.com/Shall-We-Dance/GRIP-seq

Usage:  bash GRIP-seq.sh [options]

The options include:
  --help  Print this help menu.
  --dir The root dir of our analysis
```

### Steps

1. Remove PCR deduplication

```

```

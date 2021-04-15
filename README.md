
# BVI RNA-seq Pipeline

The following instructions outline running the RNA-seq analysis pipeline for the maternal **BVI** (bacterial, virus, and immune response) project. Instructions are intended for those involved with the project at Washington University School of Medicine. As such, instructions will be provided with the assumption processing will take place on WashU RIS `compute1` compute server using `storage1` volumes. 

# Prerequisites

The following prerequisites are required for running the RNA-seq pipeline.

1. FASTQ (paired-read files)
2. Docker Image
3. Databases
4. Sample Metadata
5. Configuration File
6. Snakefile (Snakemake)

## 1. FASTQ (paired-read files)

The pipeline requires paired-end FASTQ files as input. You will supply the pipeline a file-of-filenames (fofn) listing the FASTQ read pairs, tab-delimited, one pair per line. Each FASTQ should have the fully qualified path to the file on disk listed. List all of the samples you wish to analyze in the `FASTQ.fofn` file, the minimum being one read-pair. The path to the `FASTQ.fofn` will be included in the pipeline `config.yaml` configuration file.

EXAMPLE:

```plaintext
/path/sample1.R1.fastq  /path/sample1.R2.fastq
/path/sample2.R1.fastq  /path/sample2.R2.fastq
```

## 2. Docker Image

We will be using a predefined Docker image that contains all of the required software for the pipeline. The Docker image is available through DockerHub.

[twylie/bvi_rnaseq](https://hub.docker.com/r/twylie/bvi_rnaseq)

Currently, the Docker image contains:

+ FastQC
+ Kallisto
+ MultiQC
+ Python 3.7

Note: The Docker image does not contain the required ancillary databases to run the pipeline. The database files are too large to include in the Docker image and will be hosted on `storage1` volumes.

## Databases

------------

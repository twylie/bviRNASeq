<!-- BVI RNA-seq Tutorial -->
<!-- T.N. Wylie  <twylie@wustl.edu> -->
<!-- Sun Dec 12 13:38:41 CST 2021 -->

# BVI RNA-Seq Pipeline Tutorial: RNA Clean-Up Experiment

## Overview

This tutorial walks through running the BVI RNA-seq pipeline on a set of 4 RNA double-cleanup experiments.There are 4 uBAMs in this experiment that will be processed. All files required for processing will be available to Wylie Lab members through this tutorial. For more in-depth, general details about running the BVI RNA-seq pipeline, see https://github.com/twylie/bviRNASeq.

## Collecting Sequence Data for Processing

Before we can run the RNA-seq pipeline we must locate the sequencing data and associated metadata. The information for this experimental set was provided by K.M.W. via Slack:

```plaintext
Here's the path: "/storage1/fs1/PTB/Active/2021_11_23_Test_RNA_Libraries_Remove_Small_Stuff"
Jane named the samples very clearly, and the names are in the Samplemap.csv file in that directory.
example: WLAB-RNA_1-RNA_1 had the regular clean up and WLAB-RNA_1_0-7X-RNA_1_0-7X is the same sample with a second 0.7x clean up added.
Likewise: WLAB-RNA_3-RNA_3 regular clean up and WLAB-RNA_3_0-7X-RNA_3_0-7X had the extra clean up
```

Therefore, I initially collected the uBAM (unaligned sequencing reads) files and `Samplemap.csv` (sample metadata information).

I do not like working directly with original BAM or sample info files---I usually (temporarily) copy the files I will be working with to another directory prior to processing the information. For a small set---like the 4 double clean-up experiment BAMs---this is feasible; for many (i.e. read hundreds) samples, this may not be possible, and I then recommend symbolically linking the files to another directly.

## Sequencing Read Counts and Associated Metadata

We will need to synchronize the `Samplemap.csv` and `QC.csv` metadata information with our sequencing read files in steps below, but I also like to take a precursory look beforehand to familiarize myself with the data. This information links the original uBAM, or in some cases FASTQ, files with the associated sample names. It also tells how many sequencing reads to expect per sample. I like to review the uBAMs for read count and compare to the `Samplemap.csv` information before processing. This step is not required but is a good practice.

The `Samplemap.csv` file has the following information for the tutorial samples.

| File Name                                | Flow Cell ID | Index Sequence        | Lane Number | Read Number | Sample Name                | Species | Tissue Type | Library Type | Library Name               | Completion Date | Total Bases Kb   | % PF Clusters | Avg QScore | % > Q30 | % Phix Error Rate | Protocol           |
|------------------------------------------|--------------|-----------------------|-------------|-------------|----------------------------|---------|-------------|--------------|----------------------------|-----------------|------------------|---------------|------------|---------|-------------------|--------------------|
| HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001.bam | HHYYVDSX2    | CGCTCATTAT-ACTATAGCCT |           1 | N/A         | WLAB-RNA_1-RNA_1           | human   |             | cdna library | WLAB-RNA_1-RNA_1           |      2021-11-14 | "6,719,779.652"  |         78.47 |      34.99 |    89.5 |              0.57 | Project Management |
| HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001.bam | HTMWVDSX2    | CGCTCATTAT-ACTATAGCCT |           1 | N/A         | WLAB-RNA_1_0-7X-RNA_1_0-7X | human   |             | cdna library | WLAB-RNA_1_0-7X-RNA_1_0-7X |      2021-11-18 | "26,453,648.928" |         77.76 |      34.55 |      87 |              0.68 | Project Management |
| HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001.bam | HHYYVDSX2    | GAGATTCCAT-ACTATAGCCT |           1 | N/A         | WLAB-RNA_3-RNA_3           | human   |             | cdna library | WLAB-RNA_3-RNA_3           |      2021-11-14 | "9,870,639.574"  |         78.47 |      35.73 |    92.8 |              0.57 | Project Management |
| HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004.bam | HTLTYDSX2    | GAGATTCCAT-ACTATAGCCT |           4 | N/A         | WLAB-RNA_3_0-7X-RNA_3_0-7X | human   |             | cdna library | WLAB-RNA_3_0-7X-RNA_3_0-7X |      2021-11-19 | "14,740,562.922" |         79.14 |      35.56 |    91.9 |              0.42 | Project Management |

The `QC.csv` file has the following information.

| Library                    | GTAC Sequence ID | MGI Work Order # | Date Completed | Read | Total Bases Kb (PF) | Avg QScore | % > Q30 | Phix Error Rate |
|----------------------------|------------------|------------------|----------------|------|---------------------|------------|---------|-----------------|
| WLAB-RNA_3_0-7X-RNA_3_0-7X |                  |          2865539 |                |    1 | "14,740,562.922"    |      35.86 |    93.5 |            0.27 |
| WLAB-RNA_3_0-7X-RNA_3_0-7X |                  |          2865539 |                |    2 | "14,740,562.922"    |      35.56 |    91.9 |            0.42 |
| WLAB-RNA_1-RNA_1           |                  |          2865539 |                |    1 | "6,719,779.652"     |      34.92 |    89.1 |            0.22 |
| WLAB-RNA_1-RNA_1           |                  |          2865539 |                |    2 | "6,719,779.652"     |      34.99 |    89.5 |            0.57 |
| WLAB-RNA_3-RNA_3           |                  |          2865539 |                |    1 | "9,870,639.574"     |      35.66 |    92.5 |            0.22 |
| WLAB-RNA_3-RNA_3           |                  |          2865539 |                |    2 | "9,870,639.574"     |      35.73 |    92.8 |            0.57 |
| WLAB-RNA_1_0-7X-RNA_1_0-7X |                  |          2865539 |                |    1 | "26,453,648.928"    |      35.11 |      90 |            0.33 |
| WLAB-RNA_1_0-7X-RNA_1_0-7X |                  |          2865539 |                |    2 | "26,453,648.928"    |      34.55 |      87 |            0.68 |

Right away for these data, we can see the read counts are disparate. This won’t be an issue for expression, because those values are normalized; however, keep this in mind when we look at adapter retention numbers. We will want to compare percent-of-adapter numbers among samples, etc. Also, we see between 2–5 times the target throughput of 20 million reads per sample than were originally requested.

| BAM                                      | Sample ID                  | Read Counts      |
|------------------------------------------|----------------------------|------------------|
| HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001.bam | WLAB-RNA_1-RNA_1           | 44,501,852 reads |
| HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001.bam | WLAB-RNA_1_0-7X-RNA_1_0-7X | 95,407,852 reads |
| HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001.bam | WLAB-RNA_3-RNA_3           | 65,368,474 reads |
| HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004.bam | WLAB-RNA_3_0-7X-RNA_3_0-7X | 97,619,622 reads |

## RNA-Seq Pipeline Processing Steps (Commands)

### 1. [PREREQUISITE] Convert the Input Sequencing Files to Compressed FASTQ

The pipeline requires paired-end, compressed FASTQ files as input. You will supply the pipeline a file-of-filenames (fofn) listing the FASTQ files, one filename per line, in a subsequent step.

For the tutorial data set, I was given uBAM files. We will need to convert them before running the pipeline. I will copy-in the BAM files and convert them in the following command line steps.

> **_NOTE:_** These steps have already been performed and the final 8 compressed FASTQ files are available for you here.

```plaintext
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R1_.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R2_.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001_R1_.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001_R2_.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004_R1_.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004_R2_.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R1_.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R2_.fastq.gz
```

```console
# Copy the BAMs and convert them to FASTQ files for downstream processing.

cd /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq
mkdir fastq
cd fastq
cp -pr /storage1/fs1/PTB/Active/2021_11_23_Test_RNA_Libraries_Remove_Small_Stuff/* .

# Listing after copying the BAMs.

ls -ald *

# HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001.bam
# HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001.bam
# HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004.bam
# HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001.bam
# QC.csv
# Samplemap.csv

# Here we use a Docker session that has samtools installed.

LSF_DOCKER_VOLUMES='/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq:/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq' bsub -Is -R 'select[mem>16GB] rusage[mem=16GB]' -G compute-twylie -q general-interactive -a 'docker(twylie/viromatch:latest)' zsh

# NOTE: Pairs must have the format of _R[12]_ in file name.

# Convert BAMs to FASTQ files. Note that we want the naming convention of
# R_[12].fastq.gz for processing.

samtools fastq -N -1 HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R1_.fastq -2 HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R2_.fastq HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001.bam
samtools fastq -N -1 HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001_R1_.fastq -2 HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001_R2_.fastq HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001.bam
samtools fastq -N -1 HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004_R1_.fastq -2 HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004_R2_.fastq HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004.bam
samtools fastq -N -1 HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R1_.fastq -2 HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R2_.fastq HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001.bam

# Check the FASTQ read counts and make sure they match the orginal BAM read
# counts. You may do this by looking at the Samplemap.csv read counts for the
# original, raw read counts. Compare this to your BAM and FASTQ counts.

samtools flagstat HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001.bam

# 44501852 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 mapped (0.00% : N/A)
# 44501852 + 0 paired in sequencing
# 22250926 + 0 read1
# 22250926 + 0 read2
# 0 + 0 properly paired (0.00% : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

paste - - - - < HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R[12]_.fastq | wc -l

# 44501852

# Now we will compress the FASTQ files using gzip. This can take some time, so
# consider using batch LSF jobs if you have many FASTQ files.

gzip HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R1_.fastq
gzip HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R2_.fastq
gzip HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001_R1_.fastq
gzip HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001_R2_.fastq
gzip HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004_R1_.fastq
gzip HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004_R2_.fastq
gzip HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R1_.fastq
gzip HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R2_.fastq
```

### 2. [PREREQUISITE] Locate Or Create the Transcriptome Reference File(s)

The RNA-seq pipeline uses Kallisto to map reads to the human transcriptome. Therefore, we must have a Kallisto-indexed representation of the human transcriptome. Details for creating such a resource are outlined here:

https://github.com/twylie/bviRNASeq#3-human-transcriptome-reference

NOTE: I've already created the required transcriptome reference file and associated index files here. This reference has already been indexed for Kallisto alignments.

/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/analysisReview/transcriptome_reference/Homo_sapiens.GRCh38.cdna.all.fa
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/analysisReview/transcriptome_reference/Homo_sapiens.GRCh38.cdna.all.fa.fai
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/analysisReview/transcriptome_reference/Homo_sapiens.GRCh38.cdna.all.fa.index

### 3. [PREREQUISITE] Locate Or Create the Kraken2 Database File(s)

We require a Kraken2 database setup for the bacterial characterization portion of the pipeline. More detailed instructions for creating a Kraken2 database can be found here:

https://github.com/twylie/bviRNASeq#7-kraken2-database
https://github.com/DerrickWood/kraken2/wiki/Manual

NOTE: I've already setup a "standard" DB for this purpose here:

/storage1/fs1/kwylie/Active/KRAKEN/STANDARD

### 4. [PREREQUISITE] Locate the Required Docker Image to Run the BVI RNA-Seq Pipeline

I will be using a predefined docker image that contains all of the required software for the pipeline.

NOTE: The docker image is available through dockerhub.

[[https://hub.docker.com/r/twylie/bvi_rnaseq][twylie/bvi_rnaseq]] (https://hub.docker.com/r/twylie/bvi_rnaseq)

### 5. Choose Your Directories for Processing

We will need to choose disk space for processing our data. For this tutorial, you will need to choose two areas for processing: 1) a directory on storage1, hereto referred to as the WRITE DIRECTORY; 2) a directory on scratch1, hereto refereed to as the PROCESSING DIRECTORY. You need relatively small space on scratch1 (< 1Mb) and significantly more space on storage1 (~25 Gb, not including original input BAM/FASTQ space).

When processing the data in this tutorial, I chose the following areas for my processing needs.

WRITE DIRECTORY: /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq
PROCESSING DIRECTORY: /scratch1/fs1/twylie/cleanupRNASeq

### 6. Create the FASTQ FOFN File

The pipeline requires paired-end,compressed FASTQ files as input. We will supply the pipeline a file-of-filenames (fofn) listing the FASTQ files, one filename per line. Each FASTQ should have the fully qualified path to the file on disk listed.

#+begin_src sh
# Make the FASTQ file-of-filenames file.

cd /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq
ls /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/*fastq > reads.fofn
#+end_src

This contents of the fofn file looks like this:

/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R1.fastq
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R2.fastq
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001_R1.fastq
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001_R2.fastq
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004_R1.fastq
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004_R2.fastq
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R1.fastq
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R2.fastq

### 7. Create a Sample Key for the Pipeline

We will be providing a sample key to the pipeline that associates FASTQ file paths to canonical ids. The file should be tab-delimited and contain the following three fields:

1. FASTQ Path: The fully qualified paths to the FASTQ files, which should be gzipped.
2. Canonical ID: The canonical sample ids associated with the FASTQ files. The canonical ids will be used throughout the pipeline to identify the samples.
3.  Set ID: The associated set or batch ids for the samples.

IMPORTANT! Both reads.fofn and sample_key.tsv file names must have the format of "_R[12]" to be viable.

Example of sample_key.tsv file.

FASTQ Path	Canonical ID	Set ID
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq.gz/HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R1_.fastq.gz	WLAB-RNA_1-RNA_1	1
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq.gz/HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R2_.fastq.gz	WLAB-RNA_1-RNA_1	1
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq.gz/HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R1_.fastq.gz	WLAB-RNA_1_0-7X-RNA_1_0-7X	1
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq.gz/HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R2_.fastq.gz	WLAB-RNA_1_0-7X-RNA_1_0-7X	1
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq.gz/HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001_R1_.fastq.gz	WLAB-RNA_3-RNA_3	1
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq.gz/HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001_R2_.fastq.gz	WLAB-RNA_3-RNA_3	1
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq.gz/HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004_R1_.fastq.gz	WLAB-RNA_3_0-7X-RNA_3_0-7X	1
/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq.gz/HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004_R2_.fastq.gz	WLAB-RNA_3_0-7X-RNA_3_0-7X	1

#+begin_src sh
# Creating the sample key for the placental RNA-seq samples:

cd /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq
touch sample_key.tsv

# ...now open a text editor and fill in all of the fields, while delimiting
# with tabs. CONTENT:

# FASTQ Path      Canonical ID    Set ID
# /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R1.fastq       WLAB-RNA_1-RNA_1        1
# /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HHYYVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R2.fastq       WLAB-RNA_1-RNA_1        1
# /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R1.fastq       WLAB-RNA_1_0-7X-RNA_1_0-7X      1
# /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HTMWVDSX2_CGCTCATTAT-ACTATAGCCT_L001_R2.fastq       WLAB-RNA_1_0-7X-RNA_1_0-7X      1
# /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001_R1.fastq       WLAB-RNA_3-RNA_3        1
# /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HHYYVDSX2_GAGATTCCAT-ACTATAGCCT_L001_R2.fastq       WLAB-RNA_3-RNA_3        1
# /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004_R1.fastq       WLAB-RNA_3_0-7X-RNA_3_0-7X      1
# /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/fastq/HTLTYDSX2_GAGATTCCAT-ACTATAGCCT_L004_R2.fastq       WLAB-RNA_3_0-7X-RNA_3_0-7X      1
#+end_src

### 8. Create a Directory for Processing Results

#+begin_src sh
# Create a directory to write pipeline results.

cd /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq
mkdir results
#+end_src

### 9. Create the Configuration File for Running the Pipeline

We will be passing a small configuration file that provides ancillary information for running the pipeline. The configuration file should be formatted as YAML. The configuration file will contain the following fields:

processing directory : Directory path to where the pipeline should write results.
reads fofn           : Path to the file containing a list of FASTQ paths.
transcriptome ref    : Path to the Kallisto indexed transcriptome reference file.
multiqc title        : Title text for the MultiQC web report.
multiqc description  : Description text for the MultiQC web report.
sample key           : Path to the sample key file.
kraken db            : Path to the installed/setup Kraken2 directory.

Open a text editor and fill-out the following fields. Save as config.yaml file name when finished. Your config file should look like the following file.

processing directory: '/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/results'
reads fofn: '/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/reads.fofn'
transcriptome ref: '/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/analysisReview/transcriptome_reference/Homo_sapiens.GRCh38.cdna.all.fa.index'
multiqc title: 'RNA-seq Double-Clean-Up Experiment #1'
multiqc description: 'Review of clean-up conditions for RNA-seq libraries.'
sample key: '/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/sample_key.tsv'
kraken db: '/storage1/fs1/kwylie/Active/KRAKEN/STANDARD'

Save the config as a config.yaml file in your WRITE DIRECTORY---example:

/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/config.yaml

Fields to edit include:

+ processing directory
+ reads fofn
+ multiqc title
+ multiqc description
+ sample key

### 10. Setup the Snakemake Pipeline for Processing Samples

If all of the previous tutorial steps listed above have been performed, we are ready to process the samples. We must download the Snakemake pipeline code for processing. The Snakefile (bvi_rnaseq.smk) contains the rules for running the pipeline steps. The Snakefile is static and will require no updating for running the pipeline.

There are two modes for running the pipeline: 1) single processing; 2) parallel processing. The single processing approach will run each step of the pipeline as a single event in sequential order, one at a time. Parallel processing will split steps into individual jobs and run them in parallel across the compute1 LSF server. Final output is the same for both approaches; however, parallel processing should be faster, dependent upon how many jobs were chosen for pipeline execution. I recommend always running in parallel processing mode, so the instructions listed here for this tutorial are for this processing approach.

Running the pipeline in parallel processing mode requires a little more setup. We will be adding another configuration YAML file specific for parallel processing. We will also be adding a small submitter script that helps submit individual jobs on the LSF server, via Snakemake.

We will now be working in the PROCESSING DIRECTORY on fast scratch1 space. Example:

/scratch1/fs1/twylie/cleanupRNASeq

This is where we will execute the pipeline code for processing; however, first we must download the BVI pipeline code from GitHub.

#+begin_src sh
# Download the BVI RNA-seq code from GitHub.

cd /scratch1/fs1/twylie/cleanupRNASeq
mkdir src
cd src
git clone https://github.com/twylie/bviRNASeq.git
#+end_src

We only need a few of the downloaded files for the processing our samples. Back-out to your main directory again and copy the required pipeline files to your processing directory. We will also be making a directory for LSF processing logs generated by the pipeline, as well as making the submit_lsf.py script executable.

#+begin_src sh

# Copy required pipeline files to your processing directory. Create a lsf_logs
# directory. Make sure the lsf submitter script is executable.

cd /scratch1/fs1/twylie/cleanupRNASeq
cp src/bviRNASeq/bvi_rnaseq.smk .
cp src/bviRNASeq/submit_lsf.py .
chmod +x submit_lsf.py
cp src/bviRNASeq/example/pp.yaml .
mkdir lsf_logs
#+end_src

The bvi_rnaseq.smk and submit_lsf.py files run the BVI RNA-seq pipeline. You do not have to edit these files, just make sure you've copied them in as outlined above.

The pp.yaml file is a configuration file, and will require editing for our processing needs. You may use the pp.yaml file as a template and edit it for your processing. It is very important to have both the pp.yaml and submit_lsf.py script in the same directory when processing, as the submitter script looks for the configuration file in its own directory.

Of note, Kraken2 is a memory hog. I had to use 64 GB of RAM to get it to process our tutrial RNA-seq samples. The pp.yaml file looks like this. Edit as needed.

docker:
  image: 'twylie/bvi_rnaseq'
  volumes:
    - '/scratch1/fs1/twylie'
    - '/storage1/fs1/PTB/Active'
    - '/storage1/fs1/kwylie/Active'
lsf:
  memory: '64G'
  results dir: '/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/results'
  cores: '100'
  local cores: '1'
  compute group: 'compute-kwylie'
  queue: 'general'
  latency wait: '20'
  restart times: '5'
  lsf log dir: '/scratch1/fs1/twylie/cleanupRNASeq/lsf_logs'

Fields to edit will include:

+ volumes
+ results dir
+ lsf log dir

### 11. Launch the Pipeline (Parallel Processing)

We are ready to launch the pipeline. The Snakemake command to run the pipeline looks like this.

#+begin_src sh
snakemake \
--snakefile /scratch1/fs1/twylie/cleanupRNASeq/bvi_rnaseq.smk \
--cluster /scratch1/fs1/twylie/cleanupRNASeq/submit_lsf.py \
--configfile /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/config.yaml \
--cores 100 \
--local-cores 1 \
--restart-times 5 \
--latency-wait 20 \
-p \
--rerun-incomplete
#+end_src

We can place this command in a shell script file. Example:

/scratch1/fs1/twylie/cleanupRNASeq/cmd.pp.sh

#+begin_src sh
# Make the script executable.

chmod +x /scratch1/fs1/twylie/cleanupRNASeq/cmd.pp.sh
#+end_src

#+begin_src sh
# Run the pipeline.

cd /scratch1/fs1/twylie/cleanupRNASeq
LSF_DOCKER_VOLUMES='/storage1/fs1/kwylie/Active:/storage1/fs1/kwylie/Active /scratch1/fs1/twylie:/scratch1/fs1/twylie /storage1/fs1/PTB:/storage1/fs1/PTB' bsub -M 16G -R "select[mem>16G] rusage[mem=16G]" -G compute-kwylie -q general -e $PWD/bvi.LSF.err -o $PWD/bvi.LSF.out -a 'docker(twylie/bvi_rnaseq)' sh $PWD/cmd.pp.sh
#+end_src

### 12. Results

Once finished, review results. See https://github.com/twylie/bviRNASeq#results for details on BVI RNA-seq pipeline output.

I like to move the output over to my laptop for review at this point. To do so, I run the following command on the RIS volume side.

#+begin_src sh

# Compressing the results to copy over to a laptop for review. We are skipping
# the Kraken ".out" files, as they are very large.

cd /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq
tar cvfz results.tar.gz --exclude '*.out' results/

# On the laptop side, copy from RIS to the laptop.

scp twylie@compute1-client-1.ris.wustl.edu:/storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/results.tar.gz .
tar xvfz results.tar.gz_src
cd results
#+end_src

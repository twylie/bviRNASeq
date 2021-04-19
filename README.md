
# BVI RNA-seq Pipeline

The following instructions outline running the RNA-seq analysis pipeline for the **maternal BVI** (bacterial, virus, and immune response) project. **Instructions are intended for those involved with the project at Washington University School of Medicine.** As such, instructions will be provided with the assumption processing will take place on WashU RIS `compute1` compute server using `storage1` and `scratch1` volumes.

# Author

**Todd N. Wylie** \
Assistant Professor / Department of Pediatrics \
Washington University School of Medicine \
660 S. Euclid Avenue, Campus Box 8208 \
St. Louis, MO 63110 \
Email: twylie@wustl.edu

# Prerequisites

The following prerequisite components are required for running the RNA-seq pipeline.

1. FASTQ (paired-read files)
2. Docker Image
3. Human Transcriptome Reference
4. Sample Key
5. Configuration File
6. Snakefile (Snakemake)

## 1. FASTQ (paired-read files)

The pipeline requires paired-end FASTQ files as input. You will supply the pipeline a file-of-filenames (fofn) listing the FASTQ files, one filename per line. Each FASTQ should have the fully qualified path to the file on disk listed. List all of the samples you wish to analyze in a `reads.fofn` file, the minimum being one sample; however, you must provide both read-pairs per sample. The path to the `reads.fofn` will be included in the pipeline `config.yaml` configuration file in the `reads fofn` field.

Example file in `example/reads.fofn` directory:

```plaintext
/example/HWFCGDSXX_GAATTCGT-GGCTCTGA_S25_L001_R1_001.fastq.gz
/example/HWFCGDSXX_GAATTCGT-GGCTCTGA_S25_L001_R2_001.fastq.gz
/example/AGCGATAG-AGGCGAAG_S1_R1_001.fastq.gz
/example/AGCGATAG-AGGCGAAG_S1_R2_001.fastq.gz
```

## 2. Docker Image

We will be using a predefined Docker image that contains all of the required software for the pipeline. The Docker image is available through DockerHub.

[twylie/bvi_rnaseq](https://hub.docker.com/r/twylie/bvi_rnaseq)

Currently, the Docker image contains:

+ FastQC
+ Kallisto
+ MultiQC
+ Python3
+ Pandas (Python)
+ Snakemake
+ graphviz

Note: The Docker image does not contain the required ancillary reference databases to run the pipeline. The database files are too large to include in the Docker image and will be hosted on `storage1` volumes.

## 3. Human Transcriptome Reference

We will need a human transcriptome reference for [Kallisto](https://pachterlab.github.io/kallisto/) pseudoalignments. The human transcriptome reference files are available here:

[Ensembl Database](https://www.ensembl.org/info/data/ftp/index.html)

The Kallisto authors recommend using cDNA fasta, specifically the `*.cdna.all.fa.gz` files. Kallisto can build indices directly from the gzipped files.

Example:

[http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz](http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz)

Once downloaded, we will need to index the transcriptome FASTA using Kallisto. Be sure to use the version of Kallisto in the `twylie/bvi_rnaseq` Docker Image. See [Kallisto Manual](https://pachterlab.github.io/kallisto/manual) for details on creating Kallisto indexes.

The example directory contains an example reference that has already been indexed. This example has been truncated to a single transcript (ENST00000276925.7) for testing purposes only.

```plaintext
example/transcripts.fa
example/transcripts.fa.ndx
```

## 4. Sample Key

We will be providing a sample key to the pipeline that associates FASTQ file paths to canonical ids. The file should be tab-delimited and contain the following three fields:

1. FASTQ Path: The fully qualified paths to the FASTQ files, which may be gzipped.
2. Canonical ID: The canonical sample ids associated with the FASTQ files. The canonical ids will be used throughout the pipeline to identify the samples.
3. Set ID: The associated set or batch ids for the samples.

The `example/sample_key.tsv` file provides an example sample key for reference.

```plaintext
FASTQ Path	Canonical ID	Set ID
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-AGGCGAAG_S1_R1_001.fastq.gz	TWGQ-1276-29_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-TATAGCCT_S46_R2_001.fastq.gz	TWGQ-1502-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-AGGCGAAG_S1_R2_001.fastq.gz	TWGQ-1276-29_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-AGGCGAAG_S47_R1_001.fastq.gz	TWGQ-2217-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-ATAGAGGC_S2_R1_001.fastq.gz	TWGQ-1504-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-AGGCGAAG_S47_R2_001.fastq.gz	TWGQ-2217-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-ATAGAGGC_S2_R2_001.fastq.gz	TWGQ-1504-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-ATAGAGGC_S48_R1_001.fastq.gz	TWGQ-1295-24_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-CAGGACGT_S3_R1_001.fastq.gz	TWGQ-1321-32_e_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-ATAGAGGC_S48_R2_001.fastq.gz	TWGQ-1295-24_b_R	1
...
```

The `sample_key.tsv` file should be the superset of FASTQ/ids for a cohort---i.e. it should contain all of your FASTQ files and associated canonical ids.

## 5. Configuration File

We will be passing a small configuration file that provides ancillary information for running the pipeline. The configuration file should be formatted as [YAML](https://en.wikipedia.org/wiki/YAML). The configuration file will contain the following fields:

+ processing dir: Directory path to where the pipeline should write results.
+ reads fofn: Path to the file containing a list of FASTQ paths.
+ transcriptome ref: Path to the Kallisto indexed transcriptome reference file.
+ multiqc title: Title text for the MultiQC web report.
+ multiqc description: Description text for the MultiQC web report.
+ sample key: Path to the sample key file.

The `example/config.yaml` file provides an example sample key for reference.

```YAML
processing directory: '/processing/bvi_rnaseq'
reads fofn: '/example/reads.fofn'
transcriptome ref: '/example/transcripts.fa.ndx'
multiqc title: 'BVI RNA-seq'
multiqc description: 'Initial metric review of sets #1 and #2.'
sample key: '/example/sample_key.tsv'
```

## 6. Snakefile (Snakemake)

We are using Snakemake to run our pipeline steps. A Snakefile (`bvi_rnaseq.smk`) contains the rules for running the pipeline steps. The Snakefile is static and will require no updating for running the pipeline. Currently, the pipeline involves the following rules/steps:

+ all
+ link_fastq
+ fastqc_eval
+ run_kallisto
+ run_multiqc
+ merge_kallisto_abundances
+ merge_fastqc_adapter_metrics
+ copy_multiqc_stats

<img src="images/dag.png" width="75%" height="75%" border=0 style="border:0; text-decoration:none; outline:none">

The following example rulegraph shows the workflow for two samples. Note that FastQC works on FASTQ reads pairs separately while Kallisto maps using read-pair informtion per sample.

<img src="images/rulegraph.png" width="100%" height="100%" border=0 style="border:0; text-decoration:none; outline:none">

# Test Example

## Test Data

We may use the data provided in the `example/` directory to run a simple test using the pipeline. The following test files will be used:

```plaintext
example
├── AGCGATAG-AGGCGAAG_S1_R1_001.fastq.gz
├── AGCGATAG-AGGCGAAG_S1_R2_001.fastq.gz
├── HWFCGDSXX_GAATTCGT-GGCTCTGA_S25_L001_R1_001.fastq.gz
├── HWFCGDSXX_GAATTCGT-GGCTCTGA_S25_L001_R2_001.fastq.gz
├── config.yaml
├── reads.fofn
├── sample_key.tsv
└── transcripts.fa.ndx
```

NOTE: The files listed above are simulated data and are useful only for testing purposes.

## Test Commands

```zsh
# Download the required Docker image.
docker pull twylie/bvi_rnaseq

# Clone the GitHub repository.
git clone https://github.com/twylie/bviRNASeq.git
cd bviRNAseq/

# Run an interactive Docker container.
docker container run -it -v ${PWD}:/pwd -v ${PWD}/processing:/processing -v ${PWD}/example:/example twylie/bvi_rnaseq zsh
snakemake --configfile /example/config.yaml --snakefile /pwd/bvi_rnaseq.smk --cores 1 -p

# The output of the pipeline will be in the processing directory.
ls -ald /processing/bvi_rnaseq/*
```

The top-level report files of interest will be:

+ /processing/bvi_rnaseq/BVI-RNA-seq_multiqc_report.html
+ /processing/bvi_rnaseq/abundances.merged.tsv
+ /processing/bvi_rnaseq/adapters.merged.bin70-74.tsv
+ /processing/bvi_rnaseq/multiqc_general_stats.merged.tsv

# Running Real Data

Steps are outlined here for running real data through the pipeline. Instructions are written using WashU RIS computen services, namely `storage1`, `scratch1`, and `compute1`.

## 1. Log Into WashU RIS

This step requires both VPN and RIS compute accounts at Washington University.

First, authenticate using your VPN client (e.g. Cisco AnyConnect Client). Once VPN is running successfully, we will secure-shell into RIS remotely.

```zsh
ssh twylie@compute1-client-1.ris.wustl.edu
```

## 2. Setup the Processing Directory

Because of the cache layer on `compute1` we may experience I/O latency when reading and writng many files concurrently. Therefore, we will use the faster `scratch1` space for running and writing pipeline directive files, while writing larger output files to slower, larger `storage1` space.

Setup the `storage1` processing directory first.

```zsh
cd /storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/analysisReview/
```

This will be the directory where we will be writing our analysis output from the pipeline. This directory already includes the Kallisto indexed version of the transcriptome.

```plaintext
transcriptome_reference/Homo_sapiens.GRCh38.cdna.all.fa
transcriptome_reference/Homo_sapiens.GRCh38.cdna.all.fa.fai
transcriptome_reference/Homo_sapiens.GRCh38.cdna.all.fa.index
```

We will make a results directory in this area for running the pipeline:

```zsh
mkdir /storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/analysisReview/pipelineResults/
```

This directory will be listed in the `config.yaml` file unde the `processing directory` field. Make sure there is an adequate amount of space for processing on this disk.

## 3. Setup the Working Directory

The working directory will be the area where we launch the pipeline and keep track of its associated runtime files. This diretory will be located on `scratch1` space will require a nominal amount of space. We setup a working directory here for processing and clone the GitHub repository and copy the Snakefile.

```zsh
cd /scratch1/fs1/twylie/
mkdir bviRNAseqProcessing
cd bviRNAseqProcessing/
git clone https://github.com/twylie/bviRNASeq.git
cp bviRNASeq/bvi_rnaseq.smk .
```

## 4. Make the reads.fofn File

Next we will make the `reads.fofn` that contains a list of all of the FASTQ files we wish to analyze. This file should list the fully qualified paths to the FASTQ files. FASTQ files may be gzipped. You may make this file anyway you see fit---e.g. via text editor, command line, etc.

```zsh
ls /storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/*gz > reads.fofn
grep 'SET_2' /storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/fastq.fofn >> reads.fofn
head reads.fofn
```

The `reads.fofn` file will look something like this (truncated):

```plaintext
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-AGGCGAAG_S1_R1_001.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-AGGCGAAG_S1_R2_001.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-ATAGAGGC_S2_R1_001.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-ATAGAGGC_S2_R2_001.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-CAGGACGT_S3_R1_001.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-CAGGACGT_S3_R2_001.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-CCTATCCT_S4_R1_001.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-CCTATCCT_S4_R2_001.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-GGCTCTGA_S5_R1_001.fastq.gz
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-GGCTCTGA_S5_R2_001.fastq.gz
```

NOTE: The pipeline symbolically links the FASTQ files  into the processing directory when running.

## 5. Add the sample_key.tsv File

Next we add the `sample_key.tsv` file. This file associates the canonical sample id with each FASTQ file in the `reads.fofn` file. There are only three fields in this file: (1) FASTQ Path (2) Canonical ID (3) Set ID. The Canonical ID is the unique name that will be used for a sample throughout the pipeline. Set ID is a batch number or label. If you only have one batch, just label this value as `1` in this file.

The `sample_key.tsv` file will look something like this (truncated):

```plaintext
FASTQ Path	Canonical ID	Set ID
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-AGGCGAAG_S1_R1_001.fastq.gz	TWGQ-1276-29_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-TATAGCCT_S46_R2_001.fastq.gz	TWGQ-1502-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-AGGCGAAG_S1_R2_001.fastq.gz	TWGQ-1276-29_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-AGGCGAAG_S47_R1_001.fastq.gz	TWGQ-2217-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-ATAGAGGC_S2_R1_001.fastq.gz	TWGQ-1504-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-AGGCGAAG_S47_R2_001.fastq.gz	TWGQ-2217-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-ATAGAGGC_S2_R2_001.fastq.gz	TWGQ-1504-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-ATAGAGGC_S48_R1_001.fastq.gz	TWGQ-1295-24_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-CAGGACGT_S3_R1_001.fastq.gz	TWGQ-1321-32_e_R	1
```

Make sure this file is tab-delimited and has the field labels on the first line of the file, as shown above. This file should be a superset of the FASTQ files you wish to run in the pipeline. That is, list all possible FASTQ files and associated canonical names in this file whether you run them all or not. The `reads.fofn` may include a small subset of all of the possible sample to run, but the `sample_key.tsv` file should list all of the possible samples.

```zsh
# In this case, the sample_key.tsv provided in the example/ is a viable
# metadata file for the processing I want to accomplish.
cp bviRNASeq/example/sample_key.tsv .
```

## 6. Setup the config.yaml File

Finally, we setup the `config.yaml` file. This file directs the pipeline to all of the ancillary files required to run the pipeline. The format is simple YAML. The following fields are required:

+ processing directory
+ reads fofn
+ transcriptome ref
+ multiqc title
+ multiqc description
+ sample key

This is what my `config.yaml` looks like:

```plaintext
processing directory: '/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/analysisReview/pipelineResults'
reads fofn: '/scratch1/fs1/twylie/bviRNAseqProcessing/reads.fofn'
transcriptome ref: '/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/analysisReview/transcriptome_reference/Homo_sapiens.GRCh38.cdna.all.fa.index'
multiqc title: 'BVI RNA-seq'
multiqc description: 'Maternal BVI RNA-seq Analysis (Batch #1-2)'
sample key: '/scratch1/fs1/twylie/bviRNAseqProcessing/sample_key.tsv'
```

## 7. Launch the Pipeline

There are two modes for running the pipeline: 1) single processing; 2) parallel processing. The single processing approach will run each step of the pipeline as a single event in sequential order one at a time. Parallel processing will split steps up into individual jobs and run them in parallel across the `compute1` LSF server. Final output is the same for both approaches; however, parallel processing should be faster dependent upon how many jobs were chosen for pipeline execution.

### Single Processing

The approach for single processing is fairly straight forward. Just run Snakemake as follows making sure to explicitly point to your `config.yaml` file and the main `bvi_rnaseq.smk` Snakefile.

Before submitting the job to LSF, it's easiest to place the Snakemake command statement in its own shell script, here called `cmd.sh`. The following command is placed in the shell script.

```zsh
snakemake --snakefile /scratch1/fs1/twylie/bviRNAseqProcessing/bvi_rnaseq.smk --configfile /scratch1/fs1/twylie/bviRNAseqProcessing/config.yaml --cores 1 -p
```

Once the `cmd.sh` is setup, you may launch the LSF job. Make sure you are in the working directory on `scratch1` when launching the job, in this case:

```zsh
cd /scratch1/fs1/twylie/bviRNAseqProcessing/
```

This is where Snakemake will write the `.snakemake/` directory for keeping track of Snakemake jobs.

The most verbose part of the LSF submission is making sure that Docker will see **all** of the disks that are used in the pipeline. Make sure that all of the volumes are mapped correctly in the `bsub` command. In my case, I only need two volume mappings to access all of the data I need to run the pipeline.

+ /scratch1/fs1/twylie/bviRNAseqProcessing
+ /storage1/fs1/PTB/Active

Launch the LSF job.

```zsh
LSF_DOCKER_VOLUMES='/scratch1/fs1/twylie/bviRNAseqProcessing:/scratch1/fs1/twylie/bviRNAseqProcessing /storage1/fs1/PTB/Active:/storage1/fs1/PTB/Active' \
bsub -M 16G \
-R "select[mem>16G] rusage[mem=16G]" \
-G compute-kwylie \
-q general \
-e $PWD/bvi.LSF.err \
-o $PWD/bvi.LSF.out \
-a 'docker(twylie/bvi_rnaseq)' \
sh $PWD/cmd.sh
```

The pipeline should run under LSF given the above command. To check progress:

```zsh
bjobs
```

### Parallel Processing

Running the pipeline in parallel processing mode requires a litte more setup. We will be adding another configuration YAML file specfic for parallel processing. We will also be adding a small _submitter_ script that helps submit individual jobs on the LSF server.

The LSF submission configuration YAML file looks like this:

```YAML
docker:
  image: 'twylie/bvi_rnaseq'
  volumes:
    - '/scratch1/fs1/twylie/bviRNAseqProcessing'
    - '/storage1/fs1/PTB/Active'
lsf:
  memory: '16G'
  results dir: '/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/analysisReview/pipelineResults'
  cores: '100'
  local cores: '1'
  compute group: 'compute-kwylie'
  queue: 'general'
  latency wait: '20'
  restart times: '3'
```

This information is used for each jobscript submmited to LSF.

(TO BE CONTIUNUED...)

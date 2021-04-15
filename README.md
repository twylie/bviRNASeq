
# BVI RNA-seq Pipeline

The following instructions outline running the RNA-seq analysis pipeline for the **maternal BVI** (bacterial, virus, and immune response) project. Instructions are intended for those involved with the project at Washington University School of Medicine. As such, instructions will be provided with the assumption processing will take place on WashU RIS `compute1` compute server using `storage1` and `scratch1` volumes. 

# Prerequisites

The following prerequisite components are required for running the RNA-seq pipeline.

1. FASTQ (paired-read files)
2. Docker Image
3. Human Transcriptome Reference
4. Sample Key
5. Configuration File
6. Snakefile (Snakemake)

## 1. FASTQ (paired-read files)

The pipeline requires paired-end FASTQ files as input. You will supply the pipeline a file-of-filenames (fofn) listing the FASTQ files, one filename per line. Each FASTQ should have the fully qualified path to the file on disk listed. List all of the samples you wish to analyze in a `reads.fofn` file, the minimum being one FASTQ file. The path to the `reads.fofn` will be included in the pipeline `config.yaml` configuration file in the `reads fofn` field.

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

## Reference Databases

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

## Sample Key

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
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-CAGGACGT_S3_R2_001.fastq.gz	TWGQ-1321-32_e_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-CAGGACGT_S49_R1_001.fastq.gz	TWGQ-1316-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-CCTATCCT_S4_R1_001.fastq.gz	TWGQ-2308-36_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-CAGGACGT_S49_R2_001.fastq.gz	TWGQ-1316-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-CCTATCCT_S4_R2_001.fastq.gz	TWGQ-2308-36_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-CCTATCCT_S50_R1_001.fastq.gz	TWGQ-2248-35_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-GGCTCTGA_S5_R1_001.fastq.gz	TWGQ-2139-34_e_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-CCTATCCT_S50_R2_001.fastq.gz	TWGQ-2248-35_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-GGCTCTGA_S5_R2_001.fastq.gz	TWGQ-2139-34_e_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-GGCTCTGA_S51_R1_001.fastq.gz	TWGQ-2351-14_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-TAATCTTA_S6_R1_001.fastq.gz	TWGQ-2143-14_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-GGCTCTGA_S51_R2_001.fastq.gz	TWGQ-2351-14_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-TAATCTTA_S6_R2_001.fastq.gz	TWGQ-2143-14_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-GTACTGAC_S52_R1_001.fastq.gz	TWGQ-1398-29_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-TATAGCCT_S7_R1_001.fastq.gz	TWGQ-1566-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-GTACTGAC_S52_R2_001.fastq.gz	TWGQ-1398-29_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/AGCGATAG-TATAGCCT_S7_R2_001.fastq.gz	TWGQ-1566-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-TAATCTTA_S53_R1_001.fastq.gz	TWGQ-1321-10_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-AGGCGAAG_S8_R1_001.fastq.gz	TWGQ-2326-33_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-TAATCTTA_S53_R2_001.fastq.gz	TWGQ-1321-10_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-AGGCGAAG_S8_R2_001.fastq.gz	TWGQ-2326-33_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-TATAGCCT_S54_R1_001.fastq.gz	TWGQ-1502-16_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-ATAGAGGC_S9_R1_001.fastq.gz	TWGQ-1451-15_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAATTCGT-TATAGCCT_S54_R2_001.fastq.gz	TWGQ-1502-16_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-ATAGAGGC_S9_R2_001.fastq.gz	TWGQ-1451-15_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-AGGCGAAG_S55_R1_001.fastq.gz	TWGQ-2308-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-CAGGACGT_S10_R1_001.fastq.gz	TWGQ-2232-28_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-AGGCGAAG_S55_R2_001.fastq.gz	TWGQ-2308-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-CAGGACGT_S10_R2_001.fastq.gz	TWGQ-2232-28_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-ATAGAGGC_S56_R1_001.fastq.gz	TWGQ-1543-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-CCTATCCT_S11_R1_001.fastq.gz	TWGQ-1446-19_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-ATAGAGGC_S56_R2_001.fastq.gz	TWGQ-1543-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-CCTATCCT_S11_R2_001.fastq.gz	TWGQ-1446-19_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-CAGGACGT_S57_R1_001.fastq.gz	TWGQ-1420-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-GGCTCTGA_S12_R1_001.fastq.gz	TWGQ-2188-11_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-CAGGACGT_S57_R2_001.fastq.gz	TWGQ-1420-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-GGCTCTGA_S12_R2_001.fastq.gz	TWGQ-2188-11_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-CCTATCCT_S58_R1_001.fastq.gz	TWGQ-1375-19_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-GTACTGAC_S13_R1_001.fastq.gz	TWGQ-2202-34_e_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-CCTATCCT_S58_R2_001.fastq.gz	TWGQ-1375-19_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-GTACTGAC_S13_R2_001.fastq.gz	TWGQ-2202-34_e_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-GGCTCTGA_S59_R1_001.fastq.gz	TWGQ-2139-28_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-TAATCTTA_S14_R1_001.fastq.gz	TWGQ-2202-27_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-GGCTCTGA_S59_R2_001.fastq.gz	TWGQ-2139-28_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-TAATCTTA_S14_R2_001.fastq.gz	TWGQ-2202-27_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-GTACTGAC_S60_R1_001.fastq.gz	TWGQ-1502-36_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-TATAGCCT_S15_R1_001.fastq.gz	TWGQ-2232-18_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-GTACTGAC_S60_R2_001.fastq.gz	TWGQ-1502-36_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTACTCG-TATAGCCT_S15_R2_001.fastq.gz	TWGQ-2232-18_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-TAATCTTA_S61_R1_001.fastq.gz	TWGQ-1502-32_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-AGGCGAAG_S16_R1_001.fastq.gz	TWGQ-2336-22_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-TAATCTTA_S61_R2_001.fastq.gz	TWGQ-1502-32_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-AGGCGAAG_S16_R2_001.fastq.gz	TWGQ-2336-22_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-TATAGCCT_S62_R1_001.fastq.gz	TWGQ-2312-13_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-ATAGAGGC_S17_R1_001.fastq.gz	TWGQ-1518-20_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/GAGATTCC-TATAGCCT_S62_R2_001.fastq.gz	TWGQ-2312-13_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-ATAGAGGC_S17_R2_001.fastq.gz	TWGQ-1518-20_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-AGGCGAAG_S63_R1_001.fastq.gz	TWGQ-1321-19_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-CAGGACGT_S18_R1_001.fastq.gz	TWGQ-1566-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-AGGCGAAG_S63_R2_001.fastq.gz	TWGQ-1321-19_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-CAGGACGT_S18_R2_001.fastq.gz	TWGQ-1566-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-ATAGAGGC_S64_R1_001.fastq.gz	TWGQ-1420-29_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-CCTATCCT_S19_R1_001.fastq.gz	TWGQ-2248-11_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-ATAGAGGC_S64_R2_001.fastq.gz	TWGQ-1420-29_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-CCTATCCT_S19_R2_001.fastq.gz	TWGQ-2248-11_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-CAGGACGT_S65_R1_001.fastq.gz	TWGQ-1347-27_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-GGCTCTGA_S20_R1_001.fastq.gz	TWGQ-2256-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-CAGGACGT_S65_R2_001.fastq.gz	TWGQ-1347-27_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-GGCTCTGA_S20_R2_001.fastq.gz	TWGQ-2256-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-CCTATCCT_S66_R1_001.fastq.gz	TWGQ-2326-14_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-GTACTGAC_S21_R1_001.fastq.gz	TWGQ-1344-25_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-CCTATCCT_S66_R2_001.fastq.gz	TWGQ-2326-14_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-GTACTGAC_S21_R2_001.fastq.gz	TWGQ-1344-25_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-GGCTCTGA_S67_R1_001.fastq.gz	TWGQ-2336-18_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-TAATCTTA_S22_R1_001.fastq.gz	TWGQ-1276-14_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-GGCTCTGA_S67_R2_001.fastq.gz	TWGQ-2336-18_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-TAATCTTA_S22_R2_001.fastq.gz	TWGQ-1276-14_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-GTACTGAC_S68_R1_001.fastq.gz	TWGQ-2329-29_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-TATAGCCT_S23_R1_001.fastq.gz	TWGQ-2312-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-GTACTGAC_S68_R2_001.fastq.gz	TWGQ-2329-29_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/ATTCAGAA-TATAGCCT_S23_R2_001.fastq.gz	TWGQ-2312-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-TAATCTTA_S69_R1_001.fastq.gz	TWGQ-1262-18_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-AGGCGAAG_S24_R1_001.fastq.gz	TWGQ-2351-19_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-TAATCTTA_S69_R2_001.fastq.gz	TWGQ-1262-18_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-AGGCGAAG_S24_R2_001.fastq.gz	TWGQ-2351-19_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-TATAGCCT_S70_R1_001.fastq.gz	TWGQ-1207-18_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-ATAGAGGC_S25_R1_001.fastq.gz	TWGQ-1504-16_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TAATGCGC-TATAGCCT_S70_R2_001.fastq.gz	TWGQ-1207-18_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-ATAGAGGC_S25_R2_001.fastq.gz	TWGQ-1504-16_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-AGGCGAAG_S71_R1_001.fastq.gz	TWGQ-2143-37_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-CAGGACGT_S26_R1_001.fastq.gz	TWGQ-1458-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-AGGCGAAG_S71_R2_001.fastq.gz	TWGQ-2143-37_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-CAGGACGT_S26_R2_001.fastq.gz	TWGQ-1458-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-ATAGAGGC_S72_R1_001.fastq.gz	TWGQ-1451-23_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-CCTATCCT_S27_R1_001.fastq.gz	TWGQ-1375-11_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-ATAGAGGC_S72_R2_001.fastq.gz	TWGQ-1451-23_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-CCTATCCT_S27_R2_001.fastq.gz	TWGQ-1375-11_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-CAGGACGT_S73_R1_001.fastq.gz	TWGQ-1398-25_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-GGCTCTGA_S28_R1_001.fastq.gz	TWGQ-2139-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-CAGGACGT_S73_R2_001.fastq.gz	TWGQ-1398-25_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-GGCTCTGA_S28_R2_001.fastq.gz	TWGQ-2139-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-CCTATCCT_S74_R1_001.fastq.gz	TWGQ-2256-30_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-GTACTGAC_S29_R1_001.fastq.gz	TWGQ-2217-34_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-CCTATCCT_S74_R2_001.fastq.gz	TWGQ-2256-30_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-GTACTGAC_S29_R2_001.fastq.gz	TWGQ-2217-34_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-GGCTCTGA_S75_R1_001.fastq.gz	TWGQ-2248-30_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-TAATCTTA_S30_R1_001.fastq.gz	TWGQ-2143-18_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-GGCTCTGA_S75_R2_001.fastq.gz	TWGQ-2248-30_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-TAATCTTA_S30_R2_001.fastq.gz	TWGQ-2143-18_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-TAATCTTA_S76_R1_001.fastq.gz	TWGQ-2202-15_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-TATAGCCT_S31_R1_001.fastq.gz	TWGQ-2329-23_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-TAATCTTA_S76_R2_001.fastq.gz	TWGQ-2202-15_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGCTCATT-TATAGCCT_S31_R2_001.fastq.gz	TWGQ-2329-23_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-TATAGCCT_S77_R1_001.fastq.gz	TWGQ-1458-17_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-AGGCGAAG_S32_R1_001.fastq.gz	TWGQ-1344-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGCGAA-TATAGCCT_S77_R2_001.fastq.gz	TWGQ-1458-17_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-AGGCGAAG_S32_R2_001.fastq.gz	TWGQ-1344-20_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-AGGCGAAG_S78_R1_001.fastq.gz	TWGQ-2256-24_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-ATAGAGGC_S33_R1_001.fastq.gz	TWGQ-1316-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-AGGCGAAG_S78_R2_001.fastq.gz	TWGQ-2256-24_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-ATAGAGGC_S33_R2_001.fastq.gz	TWGQ-1316-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-ATAGAGGC_S79_R1_001.fastq.gz	TWGQ-1347-15_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-CAGGACGT_S34_R1_001.fastq.gz	TWGQ-1543-16_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-ATAGAGGC_S79_R2_001.fastq.gz	TWGQ-1347-15_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-CAGGACGT_S34_R2_001.fastq.gz	TWGQ-1543-16_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-CAGGACGT_S80_R1_001.fastq.gz	TWGQ-2232-22_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-CCTATCCT_S35_R1_001.fastq.gz	TWGQ-2256-14_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-CAGGACGT_S80_R2_001.fastq.gz	TWGQ-2232-22_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-CCTATCCT_S35_R2_001.fastq.gz	TWGQ-2256-14_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-CCTATCCT_S81_R1_001.fastq.gz	TWGQ-1518-25_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-GGCTCTGA_S36_R1_001.fastq.gz	TWGQ-2188-17_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-CCTATCCT_S81_R2_001.fastq.gz	TWGQ-1518-25_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-GGCTCTGA_S36_R2_001.fastq.gz	TWGQ-2188-17_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-GGCTCTGA_S82_R1_001.fastq.gz	TWGQ-2248-19_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-TAATCTTA_S37_R1_001.fastq.gz	TWGQ-1344-16_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-GGCTCTGA_S82_R2_001.fastq.gz	TWGQ-2248-19_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-TAATCTTA_S37_R2_001.fastq.gz	TWGQ-1344-16_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-GTACTGAC_S83_R1_001.fastq.gz	TWGQ-1344-29_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-TATAGCCT_S38_R1_001.fastq.gz	TWGQ-1295-20_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-GTACTGAC_S83_R2_001.fastq.gz	TWGQ-1344-29_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CGGCTATG-TATAGCCT_S38_R2_001.fastq.gz	TWGQ-1295-20_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-TAATCTTA_S84_R1_001.fastq.gz	TWGQ-1321-27_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-AGGCGAAG_S39_R1_001.fastq.gz	TWGQ-1276-24_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-TAATCTTA_S84_R2_001.fastq.gz	TWGQ-1321-27_d_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-AGGCGAAG_S39_R2_001.fastq.gz	TWGQ-1276-24_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-TATAGCCT_S85_R1_001.fastq.gz	TWGQ-2329-18_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-ATAGAGGC_S40_R1_001.fastq.gz	TWGQ-1458-29_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCCGGAGA-TATAGCCT_S85_R2_001.fastq.gz	TWGQ-2329-18_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-ATAGAGGC_S40_R2_001.fastq.gz	TWGQ-1458-29_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCTCGCGC-ATAGAGGC_S86_R1_001.fastq.gz	TWGQ-1347-23_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-CAGGACGT_S41_R1_001.fastq.gz	TWGQ-1451-28_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCTCGCGC-ATAGAGGC_S86_R2_001.fastq.gz	TWGQ-1347-23_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-CAGGACGT_S41_R2_001.fastq.gz	TWGQ-1451-28_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCTCGCGC-CAGGACGT_S87_R1_001.fastq.gz	TWGQ-1398-20_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-CCTATCCT_S42_R1_001.fastq.gz	TWGQ-2139-23_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCTCGCGC-CAGGACGT_S87_R2_001.fastq.gz	TWGQ-1398-20_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-CCTATCCT_S42_R2_001.fastq.gz	TWGQ-2139-23_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCTCGCGC-CCTATCCT_S88_R1_001.fastq.gz	TWGQ-2351-36_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-GGCTCTGA_S43_R1_001.fastq.gz	TWGQ-2308-13_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCTCGCGC-CCTATCCT_S88_R2_001.fastq.gz	TWGQ-2351-36_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-GGCTCTGA_S43_R2_001.fastq.gz	TWGQ-2308-13_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCTCGCGC-GGCTCTGA_S89_R1_001.fastq.gz	TWGQ-2139-18_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-GTACTGAC_S44_R1_001.fastq.gz	TWGQ-2143-24_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCTCGCGC-GGCTCTGA_S89_R2_001.fastq.gz	TWGQ-2139-18_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-GTACTGAC_S44_R2_001.fastq.gz	TWGQ-2143-24_c_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCTCGCGC-TAATCTTA_S90_R1_001.fastq.gz	TWGQ-2217-14_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-TAATCTTA_S45_R1_001.fastq.gz	TWGQ-2202-10_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCTCGCGC-TAATCTTA_S90_R2_001.fastq.gz	TWGQ-2217-14_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-TAATCTTA_S45_R2_001.fastq.gz	TWGQ-2202-10_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCTCGCGC-TATAGCCT_S91_R1_001.fastq.gz	TWGQ-1420-21_b_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/CTGAAGCT-TATAGCCT_S46_R1_001.fastq.gz	TWGQ-1502-12_a_R	1
/storage1/fs1/PTB/Active/twylieAnalysis/bviRNASeq/seq/set1/TCTCGCGC-TATAGCCT_S91_R2_001.fastq.gz	TWGQ-1420-21_b_R	1
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_AGCGATAG-ATAGAGGC_S1_L001_R1_001.fastq.gz	TWGQ-2225-28_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_AGCGATAG-ATAGAGGC_S1_L001_R2_001.fastq.gz	TWGQ-2225-28_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_AGCGATAG-CCTATCCT_S2_L001_R1_001.fastq.gz	TWGQ-1299-35_d_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_AGCGATAG-CCTATCCT_S2_L001_R2_001.fastq.gz	TWGQ-1299-35_d_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_AGCGATAG-GGCTCTGA_S3_L001_R1_001.fastq.gz	TWGQ-2135-25_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_AGCGATAG-GGCTCTGA_S3_L001_R2_001.fastq.gz	TWGQ-2135-25_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_AGCGATAG-TATAGCCT_S4_L001_R1_001.fastq.gz	TWGQ-2203-12_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_AGCGATAG-TATAGCCT_S4_L001_R2_001.fastq.gz	TWGQ-2203-12_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTACTCG-ATAGAGGC_S5_L001_R1_001.fastq.gz	TWGQ-2203-16_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTACTCG-ATAGAGGC_S5_L001_R2_001.fastq.gz	TWGQ-2203-16_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTACTCG-CCTATCCT_S6_L001_R1_001.fastq.gz	TWGQ-2225-34_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTACTCG-CCTATCCT_S6_L001_R2_001.fastq.gz	TWGQ-2225-34_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTACTCG-GGCTCTGA_S7_L001_R1_001.fastq.gz	TWGQ-1311-12_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTACTCG-GGCTCTGA_S7_L001_R2_001.fastq.gz	TWGQ-1311-12_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTACTCG-TATAGCCT_S8_L001_R1_001.fastq.gz	TWGQ-2200-15_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTACTCG-TATAGCCT_S8_L001_R2_001.fastq.gz	TWGQ-2200-15_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTCAGAA-ATAGAGGC_S9_L001_R1_001.fastq.gz	TWGQ-2166-28_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTCAGAA-ATAGAGGC_S9_L001_R2_001.fastq.gz	TWGQ-2166-28_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTCAGAA-GGCTCTGA_S10_L001_R1_001.fastq.gz	TWGQ-1282-18_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTCAGAA-GGCTCTGA_S10_L001_R2_001.fastq.gz	TWGQ-1282-18_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTCAGAA-TATAGCCT_S11_L001_R1_001.fastq.gz	TWGQ-2187-15_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_ATTCAGAA-TATAGCCT_S11_L001_R2_001.fastq.gz	TWGQ-2187-15_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGCTCATT-ATAGAGGC_S12_L001_R1_001.fastq.gz	TWGQ-2203-24_d_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGCTCATT-ATAGAGGC_S12_L001_R2_001.fastq.gz	TWGQ-2203-24_d_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGCTCATT-GGCTCTGA_S13_L001_R1_001.fastq.gz	TWGQ-1311-21_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGCTCATT-GGCTCTGA_S13_L001_R2_001.fastq.gz	TWGQ-1311-21_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGCTCATT-TATAGCCT_S14_L001_R1_001.fastq.gz	TWGQ-2200-24_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGCTCATT-TATAGCCT_S14_L001_R2_001.fastq.gz	TWGQ-2200-24_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGGCTATG-ATAGAGGC_S15_L001_R1_001.fastq.gz	TWGQ-1328-27_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGGCTATG-ATAGAGGC_S15_L001_R2_001.fastq.gz	TWGQ-1328-27_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGGCTATG-CCTATCCT_S16_L001_R1_001.fastq.gz	TWGQ-1299-11_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGGCTATG-CCTATCCT_S16_L001_R2_001.fastq.gz	TWGQ-1299-11_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGGCTATG-GGCTCTGA_S17_L001_R1_001.fastq.gz	TWGQ-1366-22_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGGCTATG-GGCTCTGA_S17_L001_R2_001.fastq.gz	TWGQ-1366-22_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGGCTATG-TATAGCCT_S18_L001_R1_001.fastq.gz	TWGQ-2175-16_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CGGCTATG-TATAGCCT_S18_L001_R2_001.fastq.gz	TWGQ-2175-16_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CTGAAGCT-ATAGAGGC_S19_L001_R1_001.fastq.gz	TWGQ-1275-29_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CTGAAGCT-ATAGAGGC_S19_L001_R2_001.fastq.gz	TWGQ-1275-29_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CTGAAGCT-CCTATCCT_S20_L001_R1_001.fastq.gz	TWGQ-2135-10_a_R-2	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CTGAAGCT-CCTATCCT_S20_L001_R2_001.fastq.gz	TWGQ-2135-10_a_R-2	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CTGAAGCT-GGCTCTGA_S21_L001_R1_001.fastq.gz	TWGQ-1273-20_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CTGAAGCT-GGCTCTGA_S21_L001_R2_001.fastq.gz	TWGQ-1273-20_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CTGAAGCT-TATAGCCT_S22_L001_R1_001.fastq.gz	TWGQ-2187-23_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_CTGAAGCT-TATAGCCT_S22_L001_R2_001.fastq.gz	TWGQ-2187-23_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAATTCGT-ATAGAGGC_S23_L001_R1_001.fastq.gz	TWGQ-1275-20_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAATTCGT-ATAGAGGC_S23_L001_R2_001.fastq.gz	TWGQ-1275-20_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAATTCGT-CCTATCCT_S24_L001_R1_001.fastq.gz	TWGQ-2135-10_a_R-1	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAATTCGT-CCTATCCT_S24_L001_R2_001.fastq.gz	TWGQ-2135-10_a_R-1	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAATTCGT-GGCTCTGA_S25_L001_R1_001.fastq.gz	TWGQ-1351-19_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAATTCGT-GGCTCTGA_S25_L001_R2_001.fastq.gz	TWGQ-1351-19_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAATTCGT-TATAGCCT_S26_L001_R1_001.fastq.gz	TWGQ-2187-19_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAATTCGT-TATAGCCT_S26_L001_R2_001.fastq.gz	TWGQ-2187-19_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAGATTCC-ATAGAGGC_S27_L001_R1_001.fastq.gz	TWGQ-2166-20_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAGATTCC-ATAGAGGC_S27_L001_R2_001.fastq.gz	TWGQ-2166-20_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAGATTCC-GGCTCTGA_S28_L001_R1_001.fastq.gz	TWGQ-1311-31_d_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAGATTCC-GGCTCTGA_S28_L001_R2_001.fastq.gz	TWGQ-1311-31_d_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAGATTCC-TATAGCCT_S29_L001_R1_001.fastq.gz	TWGQ-2200-32_d_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_GAGATTCC-TATAGCCT_S29_L001_R2_001.fastq.gz	TWGQ-2200-32_d_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TAATGCGC-ATAGAGGC_S30_L001_R1_001.fastq.gz	TWGQ-1328-23_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TAATGCGC-ATAGAGGC_S30_L001_R2_001.fastq.gz	TWGQ-1328-23_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TAATGCGC-CCTATCCT_S31_L001_R1_001.fastq.gz	TWGQ-2135-10_a_R-3	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TAATGCGC-CCTATCCT_S31_L001_R2_001.fastq.gz	TWGQ-2135-10_a_R-3	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TAATGCGC-GGCTCTGA_S32_L001_R1_001.fastq.gz	TWGQ-1273-28_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TAATGCGC-GGCTCTGA_S32_L001_R2_001.fastq.gz	TWGQ-1273-28_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TAATGCGC-TATAGCCT_S33_L001_R1_001.fastq.gz	TWGQ-2187-27_d_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TAATGCGC-TATAGCCT_S33_L001_R2_001.fastq.gz	TWGQ-2187-27_d_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGCGAA-ATAGAGGC_S34_L001_R1_001.fastq.gz	TWGQ-1328-35_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGCGAA-ATAGAGGC_S34_L001_R2_001.fastq.gz	TWGQ-1328-35_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGCGAA-CCTATCCT_S35_L001_R1_001.fastq.gz	TWGQ-1299-18_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGCGAA-CCTATCCT_S35_L001_R2_001.fastq.gz	TWGQ-1299-18_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGCGAA-GGCTCTGA_S36_L001_R1_001.fastq.gz	TWGQ-1366-29_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGCGAA-GGCTCTGA_S36_L001_R2_001.fastq.gz	TWGQ-1366-29_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGCGAA-TATAGCCT_S37_L001_R1_001.fastq.gz	TWGQ-2175-20_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGCGAA-TATAGCCT_S37_L001_R2_001.fastq.gz	TWGQ-2175-20_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGGAGA-ATAGAGGC_S38_L001_R1_001.fastq.gz	TWGQ-2203-20_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGGAGA-ATAGAGGC_S38_L001_R2_001.fastq.gz	TWGQ-2203-20_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGGAGA-GGCTCTGA_S39_L001_R1_001.fastq.gz	TWGQ-1311-17_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGGAGA-GGCTCTGA_S39_L001_R2_001.fastq.gz	TWGQ-1311-17_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGGAGA-TATAGCCT_S40_L001_R1_001.fastq.gz	TWGQ-2200-20_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCCGGAGA-TATAGCCT_S40_L001_R2_001.fastq.gz	TWGQ-2200-20_b_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCTCGCGC-ATAGAGGC_S41_L001_R1_001.fastq.gz	TWGQ-2225-22_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCTCGCGC-ATAGAGGC_S41_L001_R2_001.fastq.gz	TWGQ-2225-22_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCTCGCGC-CCTATCCT_S42_L001_R1_001.fastq.gz	TWGQ-1299-29_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCTCGCGC-CCTATCCT_S42_L001_R2_001.fastq.gz	TWGQ-1299-29_c_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCTCGCGC-GGCTCTGA_S43_L001_R1_001.fastq.gz	TWGQ-2175-11_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCTCGCGC-GGCTCTGA_S43_L001_R2_001.fastq.gz	TWGQ-2175-11_a_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCTCGCGC-TATAGCCT_S44_L001_R1_001.fastq.gz	TWGQ-2175-24_d_R	2
/storage1/fs1/PTB/Active/2020_02_03_BVI_RNA_SEQ_SET_2/HWFCGDSXX_TCTCGCGC-TATAGCCT_S44_L001_R2_001.fastq.gz	TWGQ-2175-24_d_R	2

```

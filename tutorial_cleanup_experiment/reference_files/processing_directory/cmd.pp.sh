snakemake --snakefile /scratch1/fs1/twylie/cleanupRNASeq/bvi_rnaseq.smk --cluster /scratch1/fs1/twylie/cleanupRNASeq/submit_lsf.py --configfile /storage1/fs1/PTB/Active/twylieAnalysis/cleanupRNASeq/config.yaml --cores 100 --local-cores 1 --restart-times 5 --latency-wait 20 -p --rerun-incomplete
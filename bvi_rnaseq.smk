# BVI RNA-seq / Snakefile
# T.N. Wylie  <twylie@wustl.edu>
# Wed Apr 14 16:59:16 CDT 2021

import os
import re

configfile: '/example/config.yaml'


###############################################################################
#                               GLOBAL VARIABLES                              #
###############################################################################
            
OUTDIR = config['processing directory']
READS_FOFN = config['reads fofn']
SAMPLE_KEY = config['sample key']

# Get a set of the individual FASTQ files.

FASTQ = set()

with open(READS_FOFN, 'r') as fh:
    for line in fh:
        line = line.strip()
        FASTQ.add(line)

# Get a set of the individual associated canonical ids and FASTQ link names.

SAMPLE_TABLE = dict()
with open(SAMPLE_KEY, 'r') as fh:
    for i, line in enumerate(fh):
        if i == 0:
            next
        else:
            line = line.strip()
            fastq_path, canonical_id, pair_id = line.split()
            SAMPLE_TABLE.update({os.path.basename(fastq_path): canonical_id})

CANONICAL_IDS = set()
for sample in FASTQ:
    if os.path.basename(sample) in SAMPLE_TABLE.keys():
        canonical_id = SAMPLE_TABLE[os.path.basename(sample)]
        CANONICAL_IDS.add(canonical_id)

FASTQ_LINK_IDS = set()
for id_ in CANONICAL_IDS:
    FASTQ_LINK_IDS.add(id_ + '_R1')
    FASTQ_LINK_IDS.add(id_ + '_R2')

# Output directories.
    
FASTQ_OUTDIR = os.path.join(OUTDIR, 'fastq')
FASTQC_OUTDIR = os.path.join(OUTDIR, 'fastqc_results')
KALLISTO_OUTDIR = os.path.join(OUTDIR, 'kallisto_results')
MULTIQC_OUTDIR = os.path.join(OUTDIR, 'multiqc_results')


###############################################################################
#                                    RULES                                    #
###############################################################################    

rule all:
    input:
        expand(os.path.join(FASTQ_OUTDIR, '{ids}.fastq.gz'), ids = FASTQ_LINK_IDS),
        expand(os.path.join(FASTQC_OUTDIR, '{ids}_fastqc.zip'), ids = FASTQ_LINK_IDS),
        expand(os.path.join(FASTQC_OUTDIR, '{ids}_fastqc.html'), ids = FASTQ_LINK_IDS),
        expand(os.path.join(FASTQC_OUTDIR, '{ids}.done'), ids = FASTQ_LINK_IDS),
        expand(os.path.join(KALLISTO_OUTDIR, '{ids}'), ids = CANONICAL_IDS),
        expand(os.path.join(KALLISTO_OUTDIR, '{ids}/abundance.tsv'), ids = CANONICAL_IDS),
        expand(os.path.join(KALLISTO_OUTDIR, '{ids}/abundance.tsv'), ids = CANONICAL_IDS),
        expand(os.path.join(KALLISTO_OUTDIR, '{ids}/{ids}.done'), ids = CANONICAL_IDS),
        os.path.join(MULTIQC_OUTDIR, 'multiqc.done'),
        os.path.join(MULTIQC_OUTDIR, 'multiqc.log'),
        os.path.join(OUTDIR, 'abundances.merged.tsv'),
        os.path.join(OUTDIR, 'adapters.merged.bin70-74.tsv'),
        os.path.join(OUTDIR, 'multiqc_general_stats.merged.tsv')


rule link_fastq:
    input:
        samplekey = SAMPLE_KEY,
        readsfofn = READS_FOFN
    output:
        expand(os.path.join(FASTQ_OUTDIR, '{ids}.fastq.gz'), ids = FASTQ_LINK_IDS)
    log:
        os.path.join(FASTQ_OUTDIR, 'linking.log')
    run:
       cmd = ' '.join([
           'python3 /pwd/bin/link_fastq_files.py',
           '--readsfofn {input.readsfofn}',
           '--samplekey {input.samplekey}',
           '--output ' + FASTQ_OUTDIR,
           '2> {log}'
       ]) 
       shell(cmd)


rule fastqc_eval:
    input:
        os.path.join(FASTQ_OUTDIR, '{sample}.fastq.gz')
    output:
        os.path.join(FASTQC_OUTDIR, '{sample}_fastqc.zip'),
        os.path.join(FASTQC_OUTDIR, '{sample}_fastqc.html'),
        touch(os.path.join(FASTQC_OUTDIR, '{sample}.done'))
    log:
        os.path.join(FASTQC_OUTDIR, '{sample}.log'),
    run:
        cmd = ' '.join([
            'fastqc',
            '-o ' + FASTQC_OUTDIR,
            '{input}',
            '2> {log}'
        ])
        shell(cmd)


rule run_kallisto:
    input:
        r1 = os.path.join(FASTQ_OUTDIR, '{sample}_R1.fastq.gz'),
        r2 = os.path.join(FASTQ_OUTDIR, '{sample}_R2.fastq.gz')
    output:
        dir = directory(os.path.join(KALLISTO_OUTDIR, '{sample}')),
        abundance = os.path.join(KALLISTO_OUTDIR, '{sample}/abundance.tsv'),
        done = touch(os.path.join(KALLISTO_OUTDIR, '{sample}/{sample}.done'))
    log:
        os.path.join(KALLISTO_OUTDIR, '{sample}/{sample}.log')
    run:
        cmd = ' '.join([
            'kallisto quant',
            '-t 4',
            '-i ' + config['transcriptome ref'],
            '-b 50',
            '-o {output.dir}',
            '{input.r1} {input.r2}',
            '2> {log}'
        ])
        shell(cmd)


rule run_multiqc:
    input:
        expand(os.path.join(FASTQC_OUTDIR, '{ids}.done'), ids = FASTQ_LINK_IDS),
        expand(os.path.join(KALLISTO_OUTDIR, '{ids}/{ids}.done'), ids = CANONICAL_IDS)
    output:
        touch(os.path.join(MULTIQC_OUTDIR, 'multiqc.done'))
    log:
        os.path.join(MULTIQC_OUTDIR, 'multiqc.log'),
    run:
        cmd = ' '.join([
            'multiqc',
            '--title \'' + config['multiqc title'] + '\'',
            '--comment \'' + config['multiqc description'] + '\'',
            '-o ' + MULTIQC_OUTDIR,
            '-ip',
            '-v ' + ' '.join([
                KALLISTO_OUTDIR,
                FASTQC_OUTDIR
            ]),
            '2> {log}'
        ])
        shell(cmd)


rule merge_kallisto_abundances:
    input:
        expand(os.path.join(KALLISTO_OUTDIR, '{ids}/{ids}.done'), ids = CANONICAL_IDS)
    output:
        os.path.join(OUTDIR, 'abundances.merged.tsv')
    run:
        cmd = ' '.join([
            'python3 /pwd/bin/merge_kallisto_abundance.py',
            '--kallistodir ' + KALLISTO_OUTDIR,
            '--output {output}' 
        ])
        shell(cmd)


rule merge_fastqc_adapter_metrics:
    input:
        expand(os.path.join(FASTQC_OUTDIR, '{ids}.done'), ids = FASTQ_LINK_IDS)
    output:
        os.path.join(OUTDIR, 'adapters.merged.bin70-74.tsv')
    run:
        cmd = ' '.join([
            'python3 /pwd/bin/merge_fastqc_adapter.py',
            '--fastqcdir ' + FASTQC_OUTDIR,
            '--output {output}'
        ])
        shell(cmd)


rule copy_multiqc_stats:
    input:
        os.path.join(OUTDIR, 'multiqc_results/multiqc.done')
    output:
        os.path.join(OUTDIR, 'multiqc_general_stats.merged.tsv')
    run:
        cmd = ' '.join([
            'python3 /pwd/bin/copy_multiqc_stats.py',
            '--multiqc ' + os.path.join(OUTDIR, 'multiqc_results'),
            '--output {output}'
        ])
        shell(cmd)

        
# __END__

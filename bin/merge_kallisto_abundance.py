#!/usr/bin/python3.7

import argparse
import pandas as pd
import os
import re

version = '1.0'


# SUBROUTINES #################################################################

def eval_cli_arguments():

    parser = argparse.ArgumentParser(
        description='Merge multiple Kallisto abundance files.',
        prog='merge_kallisto_abundance.py',
        add_help=False
    )

    # Optional arguments.

    parser.add_argument(
        '-h',
        '--help',
        action='help',
        help='Display the extended usage statement.'
    )

    parser.add_argument(
        '--version',
        action='version',
        version=version,
        help='Display the software version number.'
    )

    # Required arguments.

    required_group = parser.add_argument_group('required')

    required_group.add_argument(
        '--kallistodir',
        metavar='DIR',
        action='store',
        help='Root path to the Kallisto results directory.',
        required=True,
    )

    required_group.add_argument(
        '--output',
        metavar='FILE',
        action='store',
        help='Path to write the merged, annotated abundance file.',
        required=True
    )

    required_group.add_argument(
        '--annotation',
        metavar='FILE',
        action='store',
        help='Path to the transcriptome annotation file.',
        required=True
    )

    return parser.parse_args()


def determine_sample_ids(arguments):

    abundance = list()

    for root, dirs, files in os.walk(arguments.kallistodir):
        for file_ in files:
            if file_.endswith('.tsv'):
                abundance.append([os.path.basename(root), '/'.join([root, file_])])

    return abundance


def merge_abundances(abundance):

    dataframes = list()

    for sample, path in abundance:
        df = pd.read_csv(path, sep='\t')
        df['sample id'] = sample
        df = df.set_index('sample id')
        dataframes.append(df)

    df_merged = pd.concat(dataframes)

    return df_merged


def write_annotated_merged_abundances(arguments):

    df_annotation = pd.read_csv(arguments.annotation, index_col='transcript_id', sep='\t')

    out = re.sub('.tsv', '.annotated.tsv', arguments.output)
    out_fh = open(out, 'w')

    out_fh.write(
        '\t'.join([
            'sample_id',
            'target_id',
            'length',
            'eff_length',
            'est_counts',
            'tpm',
            'transcript_id',
            'seq_type',
            'location',
            'gene',
            'gene_biotype',
            'transcript_biotype',
            'gene_symbol',
            'description'
        ]) + '\n'
    )

    with open(arguments.output, 'r') as fh:
        for line in fh:
            # [0] sample id
            # [1] target_id
            # [2] length
            # [3] eff_length
            # [4] est_counts
            # [5] tpm
            line = line.strip()
            abundance_fields = line.split('\t')
            target_id = abundance_fields[1]
            if target_id != 'target_id':
                # [0] transcript_id
                # [1] seq_type
                # [2] location
                # [3] gene
                # [4] gene_biotype
                # [5] transcript_biotype
                # [6] gene_symbol
                # [7] description
                annotation_fields = list(df_annotation.loc[target_id])
                annotation_fields.insert(0, target_id)
                fields = abundance_fields + annotation_fields
                out_fh.write('\t'.join(map(str, fields)) + '\n')

    return


###############################################################################
#                                     MAIN                                    #
###############################################################################

if __name__ == '__main__':

    # Given the root path to a collection of Kallisto quantification files
    # (abundance.tsv), we will create a merged table using the parent directory
    # names as the associated sample ids. Assumes that parent directory names
    # are unique.

    arguments = eval_cli_arguments()
    abundance = determine_sample_ids(arguments)
    df = merge_abundances(abundance)
    df.to_csv(arguments.output, sep='\t')
    write_annotated_merged_abundances(arguments)

# __END__

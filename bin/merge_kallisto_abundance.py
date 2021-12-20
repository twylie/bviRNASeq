#!/usr/bin/python3.7

import argparse
import pandas as pd
import os

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
        help='Path to write the merged abundance file.',
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


def annotated_abundances(arguments, df):

    df_annotation = pd.read_csv(arguments.annotation, index_col='transcript_id', sep='\t')

    annotated_abundance = dict()

    for id_ in df.index:
        target_id = df.loc[id_]['target_id']
        df_annotation.loc[target_id]  # annotation series
        df.loc[id_]  # abundance series
        abundance_dict = df.loc[id_].to_dict()
        annotation_dict = df_annotation.loc[target_id].to_dict()
        annotation_dict.update(abundance_dict)
        annotated_abundance.update({id_: annotation_dict})

    df_annotated_abundance = pd.DataFrame.from_dict(annotated_abundance).T

    return df_annotated_abundance


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
    df = annotated_abundances(arguments, df)
    df.to_csv(arguments.output, sep='\t')


# __END__

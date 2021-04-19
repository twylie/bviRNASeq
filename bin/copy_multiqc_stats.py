#!/usr/bin/python3.7

import argparse
import pandas as pd
import os
import re
import shutil

version = '1.0'


# SUBROUTINES #################################################################

def eval_cli_arguments():

    parser = argparse.ArgumentParser(
        description='Copy general stats to main directory.',
        prog='copy_multiqc_stats.py',
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
        '--multiqcdir',
        metavar='DIR',
        action='store',
        help='Root path to the MultiQC report directory.',
        required=True,
    )

    required_group.add_argument(
        '--output',
        metavar='DIR',
        action='store',
        help='Path to write the reports.',
        required=True
    )

    return parser.parse_args()


def parse_general_stats(arguments):

    # We need to locate the multiqc_general_stats.txt file path.

    for root, dirs, files in os.walk(arguments.multiqcdir):
        for file_ in files:
            if file_.endswith('multiqc_general_stats.txt'):
                stats_file = '/'.join([root, file_])
            elif file_.endswith('multiqc_report.html'):
                html_file = '/'.join([root, file_])

    # Copy the HTML report to the top-level directory.

    destination = os.path.join(os.path.dirname(arguments.output), os.path.basename(html_file))
    shutil.copyfile(html_file, destination)

    # Parse the general stats files.

    sample_ids = list()
    pair_ids = list()

    df = pd.read_csv(stats_file, sep='\t')
    for i in df.index:
        id_ = df.loc[i]['Sample']
        pairs = re.search('_R[12]', id_)
        pos = pairs.start()
        pair_id = pairs.group()[1:]
        id_ = id_[0:pos]
        sample_ids.append(id_)
        pair_ids.append(pair_id)

    df['Sample ID'] = sample_ids
    df['Pair ID'] = pair_ids

    df = df.set_index('Sample ID')

    return df


###############################################################################
#                                     MAIN                                    #
###############################################################################

if __name__ == '__main__':

    # We will copy in general stats and report files from MultiQC processing.

    arguments = eval_cli_arguments()
    df = parse_general_stats(arguments)
    df.to_csv(arguments.output, sep='\t')


# __END__

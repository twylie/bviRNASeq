import argparse
import pandas as pd
import os
import sys
import re

version = '1.0'


# SUBROUTINES #################################################################

def eval_cli_arguments():

    parser = argparse.ArgumentParser(
        description='Link FASTQ files by their canonical sample name.',
        prog='link_fastq_files.py',
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
        '--samplekey',
        metavar='FILE',
        action='store',
        help='Path to the sample key file.',
        required=True,
    )

    required_group.add_argument(
        '--readsfofn',
        metavar='FILE',
        action='store',
        help='Path to the file-of-filenames for reads.',
        required=True,
    )

    required_group.add_argument(
        '--output',
        metavar='DIR',
        action='store',
        help='Path to write FASTQ links.',
        required=True
    )

    return parser.parse_args()


def link_fastq_files(arguments, reads, df_sample_key):

    for source_path in reads:
        read_basename = os.path.basename(source_path)
        for i in df_sample_key.index:
            if re.search(read_basename, i):
                pair_id = re.search('_R[12]_', read_basename).group()
                pair_id = re.sub('_', '', pair_id)
                dest_path = os.path.join(
                    arguments.output,
                    df_sample_key.loc[i]['Canonical ID'] + '_' + pair_id + '.fastq.gz'
                )
                sys.stderr.write('[symlink] ' + source_path + ' --> ' + dest_path + '\n')
                os.symlink(source_path, dest_path)

    return


###############################################################################
#                                     MAIN                                    #
###############################################################################

if __name__ == '__main__':

    # We symbolically link in the FASTQ files by their canonical sample ids.

    arguments = eval_cli_arguments()

    reads = set()
    for read in open(arguments.readsfofn, 'r'):
        read = read.strip()
        reads.add(read)

    df_sample_key = pd.read_csv(
        arguments.samplekey,
        index_col='FASTQ Path',
        sep='\t'
    )

    link_fastq_files(arguments, reads, df_sample_key)


# __END__

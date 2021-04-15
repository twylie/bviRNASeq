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
        metavar='DIR',
        action='store',
        help='Path to write the merged abundance file.',
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


# __END__

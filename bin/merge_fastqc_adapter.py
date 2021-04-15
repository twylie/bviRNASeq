import argparse
import pandas as pd
import os
import re
import zipfile

version = '1.0'


# SUBROUTINES #################################################################

def eval_cli_arguments():

    parser = argparse.ArgumentParser(
        description='Merge multiple FastQC adapter retention reports.',
        prog='merge_fastq_adapter.py',
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
        '--fastqcdir',
        metavar='DIR',
        action='store',
        help='Root path to the FASTQC results directory.',
        required=True,
    )

    required_group.add_argument(
        '--output',
        metavar='DIR',
        action='store',
        help='Path to write the merged adapter file.',
        required=True
    )

    return parser.parse_args()


def determine_sample_ids(arguments):

    data = list()

    for root, dirs, files in os.walk(arguments.fastqcdir):
        for file_ in files:
            if file_.endswith('.zip'):
                id_ = file_
                pairs = re.search('_R[12]', id_)
                pos = pairs.start()
                pair_id = pairs.group()[1:]
                id_ = id_[0:pos]
                data.append([id_, '/'.join([root, file_]), root, pair_id])

    return data


def unzip_files(data):

    for i, [sample, path, root, pair_id] in enumerate(data):
        path_unzipped = re.sub('.zip', '', path)
        data[i].append(path_unzipped)
        with zipfile.ZipFile(path, 'r') as fh_zip:
            fh_zip.extractall(root)

    return data


def parse_fastqc_data_file(data_file_path):

    collect = False

    with open(data_file_path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if re.match('^>>Adapter Content', line):
                collect = True  # starts collection block
            if re.match('^>>END_MODULE', line):
                collect = False  # ends collection block
            if collect is True and re.match('^70-74', line):
                # Position
                # Illumina Universal Adapter
                # Illumina Small RNA 3' Adapter
                # Illumina Small RNA 5' Adapter
                # Nextera Transposase Sequence
                # SOLID Small RNA Adapter
                values = line.split('\t')

    return values


def merge_adapter_reports(data):

    dataframes = list()

    for sample, path, root, pair_id, path_unzipped in data:
        values = parse_fastqc_data_file(path_unzipped + '/fastqc_data.txt')
        values = [sample, pair_id] + values
        df = pd.DataFrame(
            [values],
            columns=[
                'Sample ID',
                'Pair ID',
                'Position',
                'Illumina Universal Adapter',
                'Illumina Small RNA 3\' Adapter',
                'Illumina Small RNA 5\' Adapter',
                'Nextera Transposase Sequence',
                'SOLID Small RNA Adapter'
            ]
        )
        df = df.set_index('Sample ID')
        dataframes.append(df)

    df_merged = pd.concat(dataframes)

    return df_merged


###############################################################################
#                                     MAIN                                    #
###############################################################################

if __name__ == '__main__':

    # Given the root path to a collection of FastQC files , we will create a
    # merged table using the parent directory names as the associated sample
    # ids. Assumes that parent directory names are unique.

    arguments = eval_cli_arguments()
    data = determine_sample_ids(arguments)
    data = unzip_files(data)
    df = merge_adapter_reports(data)
    df.to_csv(arguments.output, sep='\t')


# __END__

#!/usr/bin/python3.7

# T.N. Wylie  <twylie@wustl.edu>
# Mon Apr 19 16:31:47 CDT 2021

# This script may be used by Snakemake to submit sub-processing to the WashU
# LSF system in parallel. We will be using Docker containers to run these
# processes. The master Snakemake job will poll to see when the sub-processes
# are finished also keep track of jobs based on the original DAG.

###############################################################################
#                                DO NOT EDIT!!!                               #
###############################################################################

import os
import sys
import yaml
from snakemake.utils import read_job_properties

with open('pp.yaml', 'r') as fh:
    config = yaml.load(fh, Loader=yaml.FullLoader)

jobscript = sys.argv[-1]
props = read_job_properties(jobscript)
job_id = props['jobid']

volumes = list()
for volume in config['docker']['volumes']:
    volume += ':' + volume
    volumes.append(volume)
volumes = str(' '.join(volumes))

lsf_error = os.path.join(config['lsf']['lsf log dir'], 'LSF.err.' + str(job_id))
lsf_out = os.path.join(config['lsf']['lsf log dir'], 'LSF.out.' + str(job_id))

cmd = ' '.join([
    'LSF_DOCKER_VOLUMES="{}"'.format(volumes),
    'bsub',
    '-M {}'.format(config['lsf']['memory']),
    '-R "select[mem>{}] rusage[mem={}]"'.format(config['lsf']['memory'], config['lsf']['memory']),
    '-G {}'.format(config['lsf']['compute group']),
    '-q {}'.format(config['lsf']['queue']),
    '-e {}'.format(lsf_error),
    '-o {}'.format(lsf_out),
    '-a \'docker({})\''.format(config['docker']['image']),
    jobscript
])

print(cmd)
os.system(cmd)

# __END__

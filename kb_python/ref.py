import os
import sys

from .constants import INDEX_FILENAME, T2G_FILENAME
from .utils import create_transcript_list, run_executable


def kallisto_index(fasta_path, index_path, k=31):
    command = ['kallisto', 'index', '-i', index_path, '-k', k, fasta_path]
    run_executable(command)
    return {'index': index_path}


def create_t2g(gtf_path, t2g_path, use_name=True, use_version=False):
    r = create_transcript_list(
        gtf_path, use_name=use_name, use_version=use_version
    )
    with open(t2g_path, 'w') as f:
        for tid in r:
            if use_name:
                f.write('{}\t{}\t{}\n'.format(tid, r[tid][0], r[tid][1]))
            else:
                f.write('{}\t{}\n'.format(tid, r[tid][0]))
    return {'t2g': t2g_path}


def ref(fasta_path, gtf_path, index_path, t2g_path, overwrite=False):
    print('Generating transcript-to-gene mapping at {}'.format(t2g_path), file=sys.stderr)
    create_t2g(gtf_path, t2g_path)
    if not os.path.exists(index_path) or overwrite:
        print('Generating kallisto index at {}'.format(index_path), file=sys.stderr)
        kallisto_index(fasta_path, index_path)
    else:
        print('Skipping kallisto index because {} already exists. Use the --overwrite flag to overwrite.'.format(index_path), file=sys.stderr)

import os

from constants import (
    INDEX_FILENAME,
    REFERENCE_FILENAME,
)
from utils import (
    run_executable, )


def kallisto_index(ref_dir, k=31):
    ref_path = os.path.join(ref_dir, REFERENCE_FILENAME)
    index_path = os.path.join(ref_dir, INDEX_FILENAME)

    if not os.path.exists(index_path):
        command = ['kallisto', 'index', '-i', index_path, '-k', k, ref_path]
        run_executable(command)

    return index_path


def kallisto_bus(index_path, out_dir, technology, fastqs, threads=8):
    command = ['kallisto', 'bus']
    command += ['-i', index_path]
    command += ['-o', out_dir]
    command += ['-x', technology]
    command += ['-t', threads]
    command += fastqs
    run_executable()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', type=str, required=True)
    parser.add_argument('-x', type=str, required=True)
    parser.add_argument('-o', type=str, required=True)
    parser.add_argument('-t', type=int, default=8)
    parser.add_argument('-m', type=int, default=4)
    parser.add_argument('--h5ad', action='store_true')
    parser.add_argument('--loom', action='store_true')
    parser.add_argument('fastqs', nargs='+')
    args = parser.parse_args()

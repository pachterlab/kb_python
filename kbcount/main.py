import os

from constants import (
    BUS_FILENAME,
    BUS_S_FILENAME,
    BUS_SC_FILENAME,
    BUS_SCS_FILENAME,
    COUNTS_DIR,
    COUNTS_FILENAME,
    GENEMAP_FILENAME,
    ECMAP_FILENAME,
    INDEX_FILENAME,
    REFERENCE_FILENAME,
    TXNAMES_FILENAME,
)
from utils import (
    run_executable,
)


def kallisto_index(ref_dir, k=31):
    ref_path = os.path.join(ref_dir, REFERENCE_FILENAME)
    index_path = os.path.join(ref_dir, INDEX_FILENAME)

    if not os.path.exists(index_path):
        command = ['kallisto', 'index', '-i', index_path, '-k', k, ref_path]
        run_executable(command)

    return index_path


def kallisto_bus(fastqs, index_path, technology, out_dir, threads=8):
    command = ['kallisto', 'bus']
    command += ['-i', index_path]
    command += ['-o', out_dir]
    command += ['-x', technology]
    command += ['-t', threads]
    command += fastqs
    run_executable(command)
    return os.path.join(out_dir, BUS_FILENAME)


def bustools_sort(bus_path, out_path, threads=8, memory='4G'):
    command = ['bustools', 'sort']
    command += ['-o', out_path]
    command += ['-t', threads]
    command += ['-m', memory]
    command += [bus_path]
    run_executable(command)
    return out_path


def bustools_correct(bus_path, out_path, whitelist_path):
    command = ['bustools', 'correct']
    command += ['-o', out_path]
    command += ['-w', whitelist_path]
    command += [bus_path]
    run_executable(command)
    return out_path


def bustools_count(bus_path, out_path, genemap_path, ecmap_path, txnames_path):
    command = ['bustools', 'count']
    command += ['-o', out_path]
    command += ['-g', genemap_path]
    command += ['-e', ecmap_path]
    command += ['-t', txnames_path]
    command += ['--genecounts']
    command += [out_path]
    run_executable(command)
    return out_path


def precheck(ref_dir, technology, out_dir, fastqs, temp_dir='tmp'):
    # Check all required reference files are present.
    # Including whether the specified technology matches the provided whitelist.
    # If whitelist is not provided, use a pre-packaged one.

    # Check technololgy with number of fastqs.
    # Output what each file will be used as.

    # Setup directory structure.
    pass


def run(
        ref_dir,
        technology,
        out_dir,
        fastqs,
        temp_dir='tmp',
        threads=8,
        memory='4G'
):
    whitelist_path = 'WHITELIST'

    index_path = kallisto_index(ref_dir)
    bus_path = kallisto_bus(
        index_path, technology, out_dir, fastqs, threads=threads
    )

    bus_s_path = os.path.join(temp_dir, BUS_S_FILENAME)
    bustools_sort(bus_path, bus_s_path, threads=threads, memory=memory)
    bus_sc_path = os.path.join(temp_dir, BUS_SC_FILENAME)
    bustools_correct(bus_s_path, bus_sc_path, whitelist_path)
    bus_scs_path = os.path.join(out_dir, BUS_SCS_FILENAME)
    bustools_sort(bus_sc_path, bus_scs_path, threads=threads, memory=memory)

    counts_path = os.path.join(out_dir, COUNTS_DIR, COUNTS_FILENAME)
    genemap_path = os.path.join(ref_dir, GENEMAP_FILENAME)
    ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    txnames_path = os.path.join(out_dir, TXNAMES_FILENAME)
    bustools_count(
        bus_sc_path, counts_path, genemap_path, ecmap_path, txnames_path
    )


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', type=str, required=True)
    parser.add_argument('-x', type=str, required=True)
    parser.add_argument('-o', type=str, default='out')
    parser.add_argument('-t', type=int, default=8)
    parser.add_argument('-m', type=str, default='4G')
    parser.add_argument('--h5ad', action='store_true')
    parser.add_argument('--loom', action='store_true')
    parser.add_argument('fastqs', nargs='+')
    args = parser.parse_args()

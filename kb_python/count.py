import os
import shutil

from .constants import (
    BUS_FILENAME,
    BUS_S_FILENAME,
    BUS_SC_FILENAME,
    BUS_SCS_FILENAME,
    COUNTS_DIR,
    COUNTS_PREFIX,
    ECMAP_FILENAME,
    TXNAMES_FILENAME,
    WHITELIST_FILENAME,
)
from .utils import (
    copy_whitelist,
    get_supported_technologies,
    run_executable,
)


def kallisto_bus(fastqs, index_path, technology, out_dir, threads=8):
    command = ['kallisto', 'bus']
    command += ['-i', index_path]
    command += ['-o', out_dir]
    command += ['-x', technology]
    command += ['-t', threads]
    command += fastqs
    run_executable(command)
    return {
        'bus': os.path.join(out_dir, BUS_FILENAME),
        'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
        'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
    }


def bustools_sort(bus_path, out_path, temp_dir='tmp', threads=8, memory='4G'):
    command = ['bustools', 'sort']
    command += ['-o', out_path]
    command += ['-T', temp_dir]
    command += ['-t', threads]
    command += ['-m', memory]
    command += [bus_path]
    run_executable(command)
    return {'bus': out_path}


def bustools_correct(bus_path, out_path, whitelist_path):
    command = ['bustools', 'correct']
    command += ['-o', out_path]
    command += ['-w', whitelist_path]
    command += [bus_path]
    run_executable(command)
    return {'bus': out_path}


def bustools_count(bus_path, out_prefix, t2g_path, ecmap_path, txnames_path):
    command = ['bustools', 'count']
    command += ['-o', out_prefix]
    command += ['-g', t2g_path]
    command += ['-e', ecmap_path]
    command += ['-t', txnames_path]
    command += ['--genecounts']
    command += [bus_path]
    run_executable(command)
    return {
        'mtx': '{}.mtx'.format(out_prefix),
        'genes': '{}.genes.txt'.format(out_prefix),
        'barcodes': '{}.barcodes.txt'.format(out_prefix),
    }


def bustools_whitelist(bus_path, out_path):
    command = ['bustools', 'whitelist', '-o', out_path, bus_path]
    run_executable(command)
    return {'whitelist': out_path}


def count(
        index_path,
        t2g_path,
        technology,
        out_dir,
        fastqs,
        whitelist_path=None,
        temp_dir='tmp',
        threads=8,
        memory='4G',
        keep_temp=False,
        overwrite=False,
        loom=False,
        h5ad=False,
):
    # Make temporary directory.
    os.makedirs(temp_dir, exist_ok=True)

    bus_result = {
        'bus': os.path.join(out_dir, BUS_FILENAME),
        'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
        'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
    }
    if any(not os.path.exists(path)
           for name, path in bus_result.items()) or overwrite:
        bus_result = kallisto_bus(
            fastqs, index_path, technology, out_dir, threads=threads
        )
    else:
        print(
            'Skipping kallisto bus because output files already exist. Use the --overwrite flag to overwrite.'
        )

    sort_result = bustools_sort(
        bus_result['bus'],
        os.path.join(temp_dir, BUS_S_FILENAME),
        temp_dir=temp_dir,
        threads=threads,
        memory=memory
    )
    # Download/generate whitelist if not provided.
    if not whitelist_path:
        if technology.upper() in [t.upper()
                                  for t in get_supported_technologies()]:
            print(
                'Whitelist not provided. Copying pre-packaged {} whitelist to current directory'
                .format(technology.upper())
            )
            whitelist_path = copy_whitelist(technology)
        else:
            whitelist_path = bustools_whitelist(
                sort_result['bus'], os.path.join(out_dir, WHITELIST_FILENAME)
            )

    correct_result = bustools_correct(
        sort_result['bus'], os.path.join(temp_dir, BUS_SC_FILENAME),
        whitelist_path
    )
    sort2_result = bustools_sort(
        correct_result['bus'],
        os.path.join(out_dir, BUS_SCS_FILENAME),
        temp_dir=temp_dir,
        threads=threads,
        memory=memory
    )

    counts_dir = os.path.join(out_dir, COUNTS_DIR)
    os.makedirs(counts_dir, exist_ok=True)
    counts_prefix = os.path.join(counts_dir, COUNTS_PREFIX)
    count_result = bustools_count(
        sort2_result['bus'],
        counts_prefix,
        t2g_path,
        bus_result['ecmap'],
        bus_result['txnames'],
    )

    # Remove temporary directory.
    if not keep_temp:
        shutil.rmtree(temp_dir, ignore_errors=True)

    return count_result

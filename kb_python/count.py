import logging
import os

from .config import get_bustools_binary_path, get_kallisto_binary_path
from .constants import (
    ADATA_PREFIX,
    BUS_CDNA_PREFIX,
    BUS_FILENAME,
    BUS_INTRON_PREFIX,
    BUS_S_FILENAME,
    BUS_SC_FILENAME,
    BUS_SCS_FILENAME,
    COUNTS_DIR,
    COUNTS_PREFIX,
    ECMAP_FILENAME,
    INSPECT_FILENAME,
    TXNAMES_FILENAME,
    WHITELIST_FILENAME,
)
from .utils import (
    copy_whitelist,
    import_matrix_as_anndata,
    overlay_anndatas,
    run_executable,
    whitelist_provided,
)

logger = logging.getLogger(__name__)


def kallisto_bus(fastqs, index_path, technology, out_dir, threads=8):
    logger.info('Generating BUS file from {}'.format(' '.join(fastqs)))
    command = [get_kallisto_binary_path(), 'bus']
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
    logger.info('Sorting BUS file {} to {}'.format(bus_path, out_path))
    command = [get_bustools_binary_path(), 'sort']
    command += ['-o', out_path]
    command += ['-T', temp_dir]
    command += ['-t', threads]
    command += ['-m', memory]
    command += [bus_path]
    run_executable(command)
    return {'bus': out_path}


def bustools_inspect(bus_path, out_path, whitelist_path, ecmap_path):
    logger.info('Inspecting BUS file {}'.format(bus_path))
    command = [get_bustools_binary_path(), 'inspect']
    command += ['-o', out_path]
    command += ['-w', whitelist_path]
    command += ['-e', ecmap_path]
    command += [bus_path]
    run_executable(command)
    return {'inspect': out_path}


def bustools_correct(bus_path, out_path, whitelist_path):
    logger.info(
        'Correcting BUS records in {} to {} with whitelist {}'.format(
            bus_path, out_path, whitelist_path
        )
    )
    command = [get_bustools_binary_path(), 'correct']
    command += ['-o', out_path]
    command += ['-w', whitelist_path]
    command += [bus_path]
    run_executable(command)
    return {'bus': out_path}


def bustools_count(bus_path, out_prefix, t2g_path, ecmap_path, txnames_path):
    logger.info(
        'Generating count matrix {} from BUS file {}'.format(
            out_prefix, bus_path
        )
    )
    command = [get_bustools_binary_path(), 'count']
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


def bustools_capture(bus_path, out_path, t2c_path, ecmap_path, txnames_path):
    logger.info(
        'Capturing records from BUS file {} to {} with transcript-to-gene mapping {}'
        .format(bus_path, out_path, t2c_path)
    )
    command = [get_bustools_binary_path(), 'capture']
    command += ['-o', out_path]
    command += ['-c', t2c_path]
    command += ['-e', ecmap_path]
    command += ['-t', txnames_path]
    command += ['--transcripts']
    command += [bus_path]
    run_executable(command)
    return {'bus': bus_path}


def bustools_whitelist(bus_path, out_path):
    logger.info(
        'Generating whitelist {} from BUS file {}'.format(out_path, bus_path)
    )
    command = [
        get_bustools_binary_path(), 'whitelist', '-o', out_path, bus_path
    ]
    run_executable(command)
    return {'whitelist': out_path}


def copy_or_create_whitelist(technology, bus_path, out_dir):
    if whitelist_provided(technology):
        logger.info(
            'Copying pre-packaged {} whitelist to {}'.format(
                technology.upper(), out_dir
            )
        )
        return copy_whitelist(technology, out_dir)
    else:
        return bustools_whitelist(
            bus_path, os.path.join(out_dir, WHITELIST_FILENAME)
        )['whitelist']


def convert_matrix_to_loom(matrix_path, barcodes_path, genes_path, out_path):
    logger.info('Converting matrix {} to loom {}'.format(matrix_path, out_path))
    adata = import_matrix_as_anndata(matrix_path, barcodes_path, genes_path)
    adata.write_loom(out_path)
    return out_path


def convert_matrix_to_h5ad(matrix_path, barcodes_path, genes_path, out_path):
    logger.info('Converting matrix {} to h5ad {}'.format(matrix_path, out_path))
    adata = import_matrix_as_anndata(matrix_path, barcodes_path, genes_path)
    adata.write(out_path)
    return out_path


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
        overwrite=False,
        loom=False,
        h5ad=False,
):
    results = {}

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
        logger.info(
            'Skipping kallisto bus because output files already exist. Use the --overwrite flag to overwrite.'
        )
    results.update(bus_result)

    sort_result = bustools_sort(
        bus_result['bus'],
        os.path.join(temp_dir, BUS_S_FILENAME),
        temp_dir=temp_dir,
        threads=threads,
        memory=memory
    )
    if not whitelist_path:
        logger.info('Whitelist not provided')
        whitelist_path = copy_or_create_whitelist(
            technology, sort_result['bus'], out_dir
        )
        results.update({'whitelist': whitelist_path})

    inspect_result = bustools_inspect(
        sort_result['bus'], os.path.join(out_dir, INSPECT_FILENAME),
        whitelist_path, bus_result['ecmap']
    )
    results.update(inspect_result)
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
    results.update({'bus_scs': sort2_result['bus']})

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
    results.update(count_result)

    # Convert outputs.
    if loom:
        loom_path = convert_matrix_to_loom(
            count_result['mtx'], count_result['barcodes'],
            count_result['genes'],
            os.path.join(counts_dir, '{}.loom'.format(ADATA_PREFIX))
        )
        results.update({'loom': loom_path})
    if h5ad:
        h5ad_path = convert_matrix_to_h5ad(
            count_result['mtx'], count_result['barcodes'],
            count_result['genes'],
            os.path.join(counts_dir, '{}.h5ad'.format(ADATA_PREFIX))
        )
        results.update({'h5ad': h5ad_path})

    return results


def count_lamanno(
        index_path,
        t2g_path,
        cdna_t2c_path,
        intron_t2c_path,
        technology,
        out_dir,
        fastqs,
        whitelist_path=None,
        temp_dir='tmp',
        threads=8,
        memory='4G',
        overwrite=False,
        loom=False,
        h5ad=False,
):
    results = {}

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
        logger.info(
            'Skipping kallisto bus because output files already exist. Use the --overwrite flag to overwrite.'
        )
    results.update(bus_result)

    sort_result = bustools_sort(
        bus_result['bus'],
        os.path.join(temp_dir, BUS_S_FILENAME),
        temp_dir=temp_dir,
        threads=threads,
        memory=memory
    )
    if not whitelist_path:
        logger.info('Whitelist not provided')
        whitelist_path = copy_or_create_whitelist(
            technology, sort_result['bus'], out_dir
        )
        results.update({'whitelist': whitelist_path})

    inspect_result = bustools_inspect(
        sort_result['bus'], os.path.join(out_dir, INSPECT_FILENAME),
        whitelist_path, bus_result['ecmap']
    )
    results.update(inspect_result)
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
    results.update({'bus_scs': sort2_result['bus']})

    prefix_to_t2c = {
        BUS_CDNA_PREFIX: cdna_t2c_path,
        BUS_INTRON_PREFIX: intron_t2c_path,
    }
    counts_dir = os.path.join(out_dir, COUNTS_DIR)
    os.makedirs(counts_dir, exist_ok=True)
    for prefix, t2c_path in prefix_to_t2c.items():
        capture_result = bustools_capture(
            sort2_result['bus'],
            os.path.join(temp_dir, '{}.bus'.format(prefix)), t2c_path,
            bus_result['ecmap'], bus_result['txnames']
        )
        sort_result = bustools_sort(
            capture_result['bus'],
            os.path.join(out_dir, '{}.s.bus'.format(prefix)),
            temp_dir=temp_dir,
            threads=threads,
            memory=memory
        )

        if prefix not in results:
            results[prefix] = {}
        results[prefix].update({'bus_s': sort_result['bus']})

        counts_prefix = os.path.join(counts_dir, prefix)
        count_result = bustools_count(
            sort_result['bus'],
            counts_prefix,
            t2g_path,
            bus_result['ecmap'],
            bus_result['txnames'],
        )
        results[prefix].update(count_result)

    if loom or h5ad:
        adatas = {}
        for prefix in prefix_to_t2c:
            adatas[prefix] = import_matrix_as_anndata(
                results[prefix]['mtx'], results[prefix]['barcodes'],
                results[prefix]['genes']
            )
        adata = overlay_anndatas(adatas['spliced'], adatas['unspliced'])
        if loom:
            loom_path = os.path.join(counts_dir, '{}.loom'.format(ADATA_PREFIX))
            logger.info('Writing matrices to loom {}'.format(loom_path))
            adata.write_loom(loom_path)
            results.update({'loom': loom_path})
        if h5ad:
            h5ad_path = os.path.join(counts_dir, '{}.h5ad'.format(ADATA_PREFIX))
            logger.info('Writing matrices to h5ad {}'.format(h5ad_path))
            adata.write(h5ad_path)
            results.update({'h5ad': h5ad_path})
    return results

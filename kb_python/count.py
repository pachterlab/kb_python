import json
import os
import re
from urllib.parse import urlparse

import scipy.io

from .config import get_bustools_binary_path, get_kallisto_binary_path, is_dry
from .constants import (
    ABUNDANCE_FILENAME,
    ADATA_PREFIX,
    BATCH_FILENAME,
    BUS_CDNA_PREFIX,
    BUS_FILENAME,
    BUS_INTRON_PREFIX,
    BUS_MASHED_FILENAME,
    BUS_MERGED_FILENAME,
    CELLRANGER_BARCODES,
    CELLRANGER_DIR,
    CELLRANGER_GENES,
    CELLRANGER_MATRIX,
    CELLS_FILENAME,
    CORRECT_CODE,
    COUNTS_PREFIX,
    ECMAP_FILENAME,
    ECMAP_MERGED_FILENAME,
    FEATURE_NAME,
    FEATURE_PREFIX,
    FILTER_WHITELIST_FILENAME,
    FILTERED_CODE,
    FILTERED_COUNTS_DIR,
    GENE_NAME,
    GENES_FILENAME,
    INSPECT_FILENAME,
    KALLISTO_INFO_FILENAME,
    KB_INFO_FILENAME,
    PROJECT_CODE,
    REPORT_HTML_FILENAME,
    REPORT_NOTEBOOK_FILENAME,
    SORT_CODE,
    TCC_PREFIX,
    TXNAMES_FILENAME,
    UNFILTERED_CODE,
    UNFILTERED_COUNTS_DIR,
    WHITELIST_FILENAME,
)
from .logging import logger
from .report import render_report
from .utils import (
    copy_map,
    copy_whitelist,
    get_temporary_filename,
    import_matrix_as_anndata,
    import_tcc_matrix_as_anndata,
    make_directory,
    move_file,
    open_as_text,
    overlay_anndatas,
    read_t2g,
    run_executable,
    stream_file,
    sum_anndatas,
    update_filename,
    whitelist_provided,
)
from .stats import STATS
from .validate import validate_files

INSPECT_PARSER = re.compile(r'^.*?(?P<count>[0-9]+)')


@validate_files()
def kallisto_pseudo(batch_path, index_path, out_dir, threads=8):
    """Runs `kallisto pseudo`.

    :param batch_path: path to textfile containing batch definitions
    :type batch_path: str
    :param index_path: path to kallisto index
    :type index_path: str
    :param out_dir: path to output directory
    :type out_dir: str
    :param threads: number of threads to use, defaults to `8`
    :type threads: int, optional

    :return: dictionary containing output files
    :rtype: dict
    """
    logger.info(f'Using index {index_path} to generate matrices to {out_dir}')
    command = [get_kallisto_binary_path(), 'pseudo']
    command += ['--quant']
    command += ['-i', index_path]
    command += ['-o', out_dir]
    command += ['-b', batch_path]
    command += ['-t', threads]
    run_executable(command)

    return {
        'mtx': os.path.join(out_dir, ABUNDANCE_FILENAME),
        'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
        'cells': os.path.join(out_dir, CELLS_FILENAME),
        'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
        'info': os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    }


@validate_files()
def kallisto_bus(
    fastqs, index_path, technology, out_dir, threads=8, n=False, k=False
):
    """Runs `kallisto bus`.

    :param fastqs: list of FASTQ file paths
    :type fastqs: list
    :param index_path: path to kallisto index
    :type index_path: str
    :param technology: single-cell technology used
    :type technology: str
    :param out_dir: path to output directory
    :type out_dir: str
    :param threads: number of threads to use, defaults to `8`
    :type threads: int, optional
    :param n: include number of read in flag column (used when splitting indices),
              defaults to `False`
    :type n: bool, optional
    :param k: alignment is done per k-mer (used when splitting indices),
              defaults to `False`
    :type k: bool, optional

    :return: dictionary containing paths to generated files
    :rtype: dict
    """
    logger.info(
        f'Using index {index_path} to generate BUS file to {out_dir} from'
    )
    for fastq in fastqs:
        logger.info((' ' * 8) + fastq)
    command = [get_kallisto_binary_path(), 'bus']
    command += ['-i', index_path]
    command += ['-o', out_dir]
    command += ['-x', technology]
    command += ['-t', threads]
    if n:
        command += ['--num']
    if k:
        command += ['--kmer']
    command += fastqs
    run_executable(command)

    return {
        'bus': os.path.join(out_dir, BUS_FILENAME),
        'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
        'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
        'info': os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    }


def kallisto_bus_split(
    fastqs,
    index_paths,
    technology,
    out_dir,
    temp_dir='tmp',
    threads=8,
    memory='4G'
):
    """Runs `kallisto bus` with split indices.

    :param fastqs: list of FASTQ file paths or URLs
    :type fastqs: list
    :param index_paths: paths to kallisto indices
    :type index_paths: list
    :param technology: single-cell technology used
    :type technology: str
    :param out_dir: path to output directory
    :type out_dir: str
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional
    :param threads: number of threads to use, defaults to `8`
    :type threads: int, optional
    :param memory: amount of memory to use, defaults to `4G`
    :type memory: str, optional

    :return: dictionary containing paths to generated files
    :rtype: dict
    """
    logger.info(f'Generating BUS file using {len(index_paths)} indices')
    part_dirs = []
    for i, index_path in enumerate(index_paths):
        bus_part_dir = os.path.join(temp_dir, f'bus_part{i}')
        fastqs = stream_fastqs(fastqs, temp_dir=temp_dir)
        kallisto_bus(
            fastqs,
            index_path,
            technology,
            bus_part_dir,
            threads=threads,
            n=True,
            k=True
        )
        part_dirs.append(bus_part_dir)

        # Sort each part to temp, and then overwrite the original
        # output.bus
        bus_part = os.path.join(bus_part_dir, BUS_FILENAME)
        sort_part = bustools_sort(
            bus_part,
            get_temporary_filename(temp_dir),
            temp_dir=temp_dir,
            threads=threads,
            memory=memory,
            flags=True
        )
        move_file(sort_part['bus'], bus_part)

    # Mash parts into one & sort by flag again
    mash_result = bustools_mash(part_dirs, out_dir)
    sort_result = bustools_sort(
        mash_result['bus'],
        get_temporary_filename(temp_dir),
        temp_dir=temp_dir,
        threads=threads,
        memory=memory,
        flags=True
    )
    move_file(sort_result['bus'], mash_result['bus'])

    # Merge
    merge_result = bustools_merge(
        mash_result['bus'], out_dir, mash_result['ecmap'],
        mash_result['txnames']
    )

    # Move files to appropriate places
    bus_path = os.path.join(out_dir, BUS_FILENAME)
    ecmap_path = os.path.join(out_dir, ECMAP_FILENAME)
    move_file(merge_result['bus'], bus_path)
    move_file(merge_result['ecmap'], ecmap_path)

    return {
        'bus': bus_path,
        'ecmap': ecmap_path,
        'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
        'info': os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    }


@validate_files(pre=False)
def bustools_mash(out_dirs, out_dir):
    """Runs `bustools mash`. Additionally, combines the `run_info.json`s into
    one.

    :param out_dirs: list of `kallisto bus` output directories. Note that
                     BUS files should be sorted by flag
    :type out_dirs: list
    :param out_dir: path to output directory
    :type out_dir: str

    :return: dictionary containing paths to generated files
    :rtype: dict
    """
    logger.info(f'Mashing BUS records to {out_dir} from')
    for out in out_dirs:
        logger.info((' ' * 8) + out)
    command = [get_bustools_binary_path(), 'mash']
    command += ['-o', out_dir]
    command += out_dirs
    run_executable(command)

    # Combine run_info.jsons (don't run this if dry)
    info_path = None
    if not is_dry():
        run_info = {}
        for o_dir in out_dirs:
            info_path = os.path.join(o_dir, KALLISTO_INFO_FILENAME)
            with open(info_path, 'r') as f:
                info = json.load(f)

            for key, value in info.items():
                run_info.setdefault(key, []).append(value)
        info_path = os.path.join(out_dir, KALLISTO_INFO_FILENAME)
        with open(info_path, 'w') as f:
            json.dump(run_info, f, indent=4)
    return {
        'bus': os.path.join(out_dir, BUS_MASHED_FILENAME),
        'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
        'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
        'info': info_path
    }


@validate_files(pre=False)
def bustools_merge(bus_path, out_dir, ecmap_path, txnames_path):
    """Runs `bustools merge`.

    :param bus_path: path to BUS file to merge
    :type bus_path: str
    :param out_dir: path to output directory, where the merged BUS file and
                    ecmap will be written
    :type out_dir: str
    :param ecmap_path: path to ecmap file, as generated by `kallisto bus`
    :type ecmap_path: str
    :param txnames_path: path to transcript names file, as generated by `kallisto bus`
    :type txnames_path: str

    :return: dictionary containing path to generated BUS file and merged ecmap
    :rtype: dict
    """
    logger.info(f'Merging BUS records in {bus_path} to {out_dir}')
    command = [get_bustools_binary_path(), 'merge']
    command += ['-o', out_dir]
    command += ['-e', ecmap_path]
    command += ['-t', txnames_path]
    command += [bus_path]
    run_executable(command)

    return {
        'bus': os.path.join(out_dir, BUS_MERGED_FILENAME),
        'ecmap': os.path.join(out_dir, ECMAP_MERGED_FILENAME),
    }


@validate_files(pre=False)
def bustools_project(bus_path, out_path, map_path, ecmap_path, txnames_path):
    """Runs `bustools project`.

    :param bus_path: path to BUS file to sort
    :type bus_path: str
    :param out_dir: path to output directory
    :type out_dir: str
    :param map_path: path to file containing source-to-destination mapping
    :type map_path: str
    :param ecmap_path: path to ecmap file, as generated by `kallisto bus`
    :type ecmap_path: str
    :param txnames_path: path to transcript names file, as generated by `kallisto bus`
    :type txnames_path: str

    :return: dictionary containing path to generated BUS file
    :rtype: dict
    """
    logger.info('Projecting BUS file {} with map {}'.format(bus_path, map_path))
    command = [get_bustools_binary_path(), 'project']
    command += ['-o', out_path]
    command += ['-m', map_path]
    command += ['-e', ecmap_path]
    command += ['-t', txnames_path]
    command += ['--barcode']
    command += [bus_path]
    run_executable(command)
    return {'bus': out_path}


@validate_files(pre=False)
def bustools_sort(
    bus_path, out_path, temp_dir='tmp', threads=8, memory='4G', flags=False
):
    """Runs `bustools sort`.

    :param bus_path: path to BUS file to sort
    :type bus_path: str
    :param out_dir: path to output BUS path
    :type out_dir: str
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional
    :param threads: number of threads to use, defaults to `8`
    :type threads: int, optional
    :param memory: amount of memory to use, defaults to `4G`
    :type memory: str, optional
    :param flags: whether to supply the `--flags` argument to sort, defaults to
                  `False`
    :type flags: bool, optional

    :return: dictionary containing path to generated index
    :rtype: dict
    """
    logger.info('Sorting BUS file {} to {}'.format(bus_path, out_path))
    command = [get_bustools_binary_path(), 'sort']
    command += ['-o', out_path]
    command += ['-T', temp_dir]
    command += ['-t', threads]
    command += ['-m', memory]
    if flags:
        command += ['--flags']
    command += [bus_path]
    run_executable(command)
    return {'bus': out_path}


@validate_files(pre=False)
def bustools_inspect(bus_path, out_path, whitelist_path, ecmap_path):
    """Runs `bustools inspect`.

    :param bus_path: path to BUS file to sort
    :type bus_path: str
    :param out_path: path to output inspect JSON file
    :type out_path: str
    :param whitelist_path: path to whitelist
    :type whitelist_path: str
    :param ecmap_path: path to ecmap file, as generated by `kallisto bus`
    :type ecmap_path: str

    :return: dictionary containing path to generated index
    :rtype: dict
    """
    logger.info('Inspecting BUS file {}'.format(bus_path))
    command = [get_bustools_binary_path(), 'inspect']
    command += ['-o', out_path]
    command += ['-w', whitelist_path]
    command += ['-e', ecmap_path]
    command += [bus_path]
    run_executable(command)
    return {'inspect': out_path}


@validate_files(pre=False)
def bustools_correct(bus_path, out_path, whitelist_path):
    """Runs `bustools correct`.

    :param bus_path: path to BUS file to correct
    :type bus_path: str
    :param out_path: path to output corrected BUS file
    :type out_path: str
    :param whitelist_path: path to whitelist
    :type whitelist_path: str

    :return: dictionary containing path to generated index
    :rtype: dict
    """
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


@validate_files(pre=False)
def bustools_count(
    bus_path,
    out_prefix,
    t2g_path,
    ecmap_path,
    txnames_path,
    tcc=False,
    mm=False
):
    """Runs `bustools count`.

    :param bus_path: path to BUS file to correct
    :type bus_path: str
    :param out_prefix: prefix of the output files to generate
    :type out_prefix: str
    :param t2g_path: path to output transcript-to-gene mapping
    :type t2g_path: str
    :param ecmap_path: path to ecmap file, as generated by `kallisto bus`
    :type ecmap_path: str
    :param txnames_path: path to transcript names file, as generated by `kallisto bus`
    :type txnames_path: str
    :param tcc: whether to generate a TCC matrix instead of a gene count matrix,
                defaults to `False`
    :type tcc: bool, optional
    :param mm: whether to include BUS records that pseudoalign to multiple genes,
               defaults to `False`
    :type mm: bool, optional

    :return: dictionary containing path to generated index
    :rtype: dict
    """
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
    if not tcc:
        command += ['--genecounts']
    if mm:
        command += ['--multimapping']
    command += [bus_path]
    run_executable(command)
    return {
        'mtx':
            '{}.mtx'.format(out_prefix),
        'ec' if tcc else 'genes':
            '{}.ec.txt'.format(out_prefix)
            if tcc else '{}.genes.txt'.format(out_prefix),
        'barcodes':
            '{}.barcodes.txt'.format(out_prefix),
    }


@validate_files(pre=False)
def bustools_capture(
    bus_path,
    out_path,
    capture_path,
    ecmap_path,
    txnames_path,
    capture_type='transcripts'
):
    """Runs `bustools capture`.

    :param bus_path: path to BUS file to capture
    :type bus_path: str
    :param out_path: path to BUS file to generate
    :type out_path: str
    :param capture_path: path transcripts-to-capture list
    :type capture_path: str
    :param ecmap_path: path to ecmap file, as generated by `kallisto bus`
    :type ecmap_path: str
    :param txnames_path: path to transcript names file, as generated by `kallisto bus`
    :type txnames_path: str
    :param capture_type: the type of information in the capture list.
                      can be one of `transcripts`, `umis`, `barcode`.
    :type capture_type: str

    :return: dictionary containing path to generated index
    :rtype: dict
    """
    logger.info(
        'Capturing records from BUS file {} to {} with capture list {}'.format(
            bus_path, out_path, capture_path
        )
    )
    command = [get_bustools_binary_path(), 'capture']
    command += ['-o', out_path]
    command += ['-c', capture_path]
    command += ['-e', ecmap_path]
    command += ['-t', txnames_path]
    command += ['--complement']
    command += ['--{}'.format(capture_type)]
    command += [bus_path]
    run_executable(command)
    return {'bus': out_path}


@validate_files(pre=False)
def bustools_whitelist(bus_path, out_path):
    """Runs `bustools whitelist`.

    :param bus_path: path to BUS file generate the whitelist from
    :type bus_path: str
    :param out_path: path to output whitelist
    :type out_path: str

    :return: dictionary containing path to generated index
    :rtype: dict
    """
    logger.info(
        'Generating whitelist {} from BUS file {}'.format(out_path, bus_path)
    )
    command = [
        get_bustools_binary_path(), 'whitelist', '-o', out_path, bus_path
    ]
    run_executable(command)
    return {'whitelist': out_path}


def write_smartseq_batch(fastq_pairs, cell_ids, out_path):
    """Write a 3-column TSV specifying batch information for Smart-seq reads.
    This file is required to use `kallisto pseudo` on multiple samples (= cells).

    :param fastq_pairs: list of pairs of FASTQs
    :type fastq_pairs: list
    :param cell_ids: list of cell IDs
    :type cell_ids: list
    :param out_path: path to batch file to output
    :type out_path: str

    :return: dictionary of written batch file
    :rtype: dict
    """
    logger.info(f'Writing batch definition TSV to {out_path}')
    with open(out_path, 'w') as f:
        for cell_id, (fastq_1, fastq_2) in zip(cell_ids, fastq_pairs):
            f.write(f'{cell_id}\t{fastq_1}\t{fastq_2}\n')
    return {'batch': out_path}


def matrix_to_cellranger(
    matrix_path, barcodes_path, genes_path, t2g_path, out_dir
):
    """Convert bustools count matrix to cellranger-format matrix.

    :param matrix_path: path to matrix
    :type matrix_path: str
    :param barcodes_path: list of paths to barcodes.txt
    :type barcodes_path: str
    :param genes_path: path to genes.txt
    :type genes_path: str
    :param t2g_path: path to transcript-to-gene mapping
    :type t2g_path: str
    :param out_dir: path to output matrix
    :type out_dir: str

    :return: dictionary of matrix files
    :rtype: dict
    """
    make_directory(out_dir)
    logger.info(f'Writing matrix in cellranger format to {out_dir}')

    cr_matrix_path = os.path.join(out_dir, CELLRANGER_MATRIX)
    cr_barcodes_path = os.path.join(out_dir, CELLRANGER_BARCODES)
    cr_genes_path = os.path.join(out_dir, CELLRANGER_GENES)

    # Cellranger outputs genes x cells matrix
    mtx = scipy.io.mmread(matrix_path)
    scipy.io.mmwrite(cr_matrix_path, mtx.T, field='integer')

    with open(barcodes_path, 'r') as f, open(cr_barcodes_path, 'w') as out:
        for line in f:
            if line.isspace():
                continue
            out.write(f'{line.strip()}-1\n')

    # Get all (available) gene names
    gene_to_name = {}
    with open(t2g_path, 'r') as f:
        for line in f:
            if line.isspace():
                continue
            split = line.strip().split('\t')
            if len(split) > 2:
                gene_to_name[split[1]] = split[2]

    with open(genes_path, 'r') as f, open(cr_genes_path, 'w') as out:
        for line in f:
            if line.isspace():
                continue
            gene = line.strip()
            gene_name = gene_to_name.get(gene, gene)
            out.write(f'{gene}\t{gene_name}\n')

    return {
        'mtx': cr_matrix_path,
        'barcodes': cr_barcodes_path,
        'genes': cr_genes_path
    }


def convert_matrix(
    counts_dir,
    matrix_path,
    barcodes_path,
    genes_path=None,
    ec_path=None,
    t2g_path=None,
    txnames_path=None,
    name='gene',
    loom=False,
    h5ad=False,
    tcc=False,
    threads=8,
):
    """Convert a gene count or TCC matrix to loom or h5ad.

    :param counts_dir: path to counts directory
    :type counts_dir: str
    :param matrix_path: path to matrix
    :type matrix_path: str
    :param barcodes_path: list of paths to barcodes.txt
    :type barcodes_path: str
    :param genes_path: path to genes.txt, defaults to `None`
    :type genes_path: str, optional
    :param ec_path: path to ec.txt, defaults to `None`
    :type ec_path: str, optional
    :param t2g_path: path to transcript-to-gene mapping. If this is provided,
                     the third column of the mapping is appended to the
                     anndata var, defaults to `None`
    :type t2g_path: str, optional
    :param txnames_path: path to transcripts.txt, defaults to `None`
    :type txnames_path: str, optional
    :param name: name of the columns, defaults to "gene"
    :type name: str, optional
    :param loom: whether to generate loom file, defaults to `False`
    :type loom: bool, optional
    :param h5ad: whether to generate h5ad file, defaults to `False`
    :type h5ad: bool, optional
    :param tcc: whether the matrix is a TCC matrix, defaults to `False`
    :type tcc: bool, optional
    :param threads: number of threads to use, defaults to `8`
    :type threads: int, optional

    :return: dictionary of generated files
    :rtype: dict
    """
    results = {}
    logger.info('Reading matrix {}'.format(matrix_path))
    adata = import_tcc_matrix_as_anndata(
        matrix_path, barcodes_path, ec_path, txnames_path, threads=threads
    ) if tcc else import_matrix_as_anndata(
        matrix_path, barcodes_path, genes_path, t2g_path=t2g_path, name=name
    )
    if loom:
        loom_path = os.path.join(counts_dir, '{}.loom'.format(ADATA_PREFIX))
        logger.info('Writing matrix to loom {}'.format(loom_path))
        adata.write_loom(loom_path)
        results.update({'loom': loom_path})
    if h5ad:
        h5ad_path = os.path.join(counts_dir, '{}.h5ad'.format(ADATA_PREFIX))
        logger.info('Writing matrix to h5ad {}'.format(h5ad_path))
        adata.write(h5ad_path)
        results.update({'h5ad': h5ad_path})

    return results


def convert_matrices(
    counts_dir,
    matrix_paths,
    barcodes_paths,
    genes_paths=None,
    ec_paths=None,
    t2g_path=None,
    txnames_path=None,
    name='gene',
    loom=False,
    h5ad=False,
    nucleus=False,
    tcc=False,
    threads=8,
):
    """Convert a gene count or TCC matrix to loom or h5ad.

    :param counts_dir: path to counts directory
    :type counts_dir: str
    :param matrix_paths: list of paths to matrices
    :type matrix_paths: list
    :param barcodes_paths: list of paths to barcodes.txt
    :type barcodes_paths: list
    :param genes_paths: list of paths to genes.txt, defaults to `None`
    :type genes_paths: list, optional
    :param ec_paths: list of path to ec.txt, defaults to `None`
    :type ec_paths: list, optional
    :param t2g_path: path to transcript-to-gene mapping. If this is provided,
                     the third column of the mapping is appended to the
                     anndata var, defaults to `None`
    :type t2g_path: str, optional
    :param txnames_path: list of paths to transcripts.txt, defaults to `None`
    :type txnames_path: str, optional
    :param name: name of the columns, defaults to "gene"
    :type name: str, optional
    :param loom: whether to generate loom file, defaults to `False`
    :type loom: bool, optional
    :param h5ad: whether to generate h5ad file, defaults to `False`
    :type h5ad: bool, optional
    :param nucleus: whether the matrices contain single nucleus counts, defaults to `False`
    :type nucleus: bool, optional
    :param tcc: whether the matrix is a TCC matrix, defaults to `False`
    :type tcc: bool, optional
    :param threads: number of threads to use, defaults to `8`
    :type threads: int, optional

    :return: dictionary of generated files
    :rtype: dict
    """
    results = {}
    adatas = []
    matrix_paths = matrix_paths or []
    barcodes_paths = barcodes_paths or []
    genes_paths = genes_paths or []
    ec_paths = ec_paths or []
    for matrix_path, barcodes_path, genes_ec_path in zip(
            matrix_paths, barcodes_paths, ec_paths
            if not genes_paths or None in genes_paths else genes_paths):
        logger.info('Reading matrix {}'.format(matrix_path))
        adatas.append(
            import_tcc_matrix_as_anndata(
                matrix_path,
                barcodes_path,
                genes_ec_path,
                txnames_path,
                threads=threads
            ) if tcc else import_matrix_as_anndata(
                matrix_path,
                barcodes_path,
                genes_ec_path,
                t2g_path=t2g_path,
                name=name
            )
        )
    logger.info('Combining matrices')
    adata = sum_anndatas(*adatas) if nucleus else overlay_anndatas(*adatas)
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


def filter_with_bustools(
    bus_path,
    ecmap_path,
    txnames_path,
    t2g_path,
    whitelist_path,
    filtered_bus_path,
    counts_prefix=None,
    tcc=False,
    mm=False,
    kite=False,
    temp_dir='tmp',
    threads=8,
    memory='4G',
    count=True,
    loom=False,
    h5ad=False,
    cellranger=False
):
    """Generate filtered count matrices with bustools.

    :param bus_path: path to sorted, corrected, sorted BUS file
    :type bus_path: str
    :param ecmap_path: path to matrix ec file
    :type ecmap_path: str
    :param txnames_path: path to list of transcripts
    :type txnames_path: str
    :param t2g_path: path to transcript-to-gene mapping
    :type t2g_path: str
    :param whitelist_path: path to filter whitelist to generate
    :type whitelist_path: str
    :param filtered_bus_path: path to filtered BUS file to generate
    :type filtered_bus_path: str
    :param counts_prefix: prefix of count matrix, defaults to `None`
    :type counts_prefix: str, optional
    :param tcc: whether to generate a TCC matrix instead of a gene count matrix,
                defaults to `False`
    :type tcc: bool, optional
    :param mm: whether to include BUS records that pseudoalign to multiple genes,
               defaults to `False`
    :type mm: bool, optional
    :param kite: Whether this is a KITE workflow
    :type kite: bool, optional
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional
    :param threads: number of threads to use, defaults to `8`
    :type threads: int, optional
    :param memory: amount of memory to use, defaults to `4G`
    :type memory: str, optional
    :param loom: whether to convert the final count matrix into a loom file,
                 defaults to `False`
    :type loom: bool, optional
    :param h5ad: whether to convert the final count matrix into a h5ad file,
                 defaults to `False`
    :type h5ad: bool, optional
    :param cellranger: whether to convert the final count matrix into a
                       cellranger-compatible matrix, defaults to `False`
    :type cellranger: bool, optional

    :return: dictionary of generated files
    :rtype: dict
    """
    logger.info('Filtering with bustools')
    results = {}
    whitelist_result = bustools_whitelist(bus_path, whitelist_path)
    results.update(whitelist_result)
    correct_result = bustools_correct(
        bus_path,
        os.path.join(
            temp_dir, update_filename(os.path.basename(bus_path), CORRECT_CODE)
        ),
        whitelist_result['whitelist'],
    )
    sort_result = bustools_sort(
        correct_result['bus'],
        filtered_bus_path,
        temp_dir=temp_dir,
        threads=threads,
        memory=memory,
    )
    results.update({'bus_scs': sort_result['bus']})

    if count:
        counts_dir = os.path.dirname(counts_prefix)
        make_directory(counts_dir)
        count_result = bustools_count(
            sort_result['bus'],
            counts_prefix,
            t2g_path,
            ecmap_path,
            txnames_path,
            tcc=tcc,
            mm=mm,
        )
        results.update(count_result)

        if loom or h5ad:
            results.update(
                convert_matrix(
                    counts_dir,
                    count_result['mtx'],
                    count_result['barcodes'],
                    genes_path=count_result.get('genes'),
                    t2g_path=t2g_path,
                    ec_path=count_result.get('ec'),
                    txnames_path=txnames_path,
                    name=FEATURE_NAME if kite else GENE_NAME,
                    loom=loom,
                    h5ad=h5ad,
                    tcc=tcc,
                    threads=threads
                )
            )
        if cellranger:
            if not tcc:
                cr_result = matrix_to_cellranger(
                    count_result['mtx'], count_result['barcodes'],
                    count_result['genes'], t2g_path,
                    os.path.join(counts_dir, CELLRANGER_DIR)
                )
                results.update({'cellranger': cr_result})
            else:
                logger.warning(
                    'TCC matrices can not be converted to cellranger-compatible format.'
                )

    return results


def stream_fastqs(fastqs, temp_dir='tmp'):
    """Given a list of fastqs (that may be local or remote paths), stream any
    remote files. Internally, calls utils.

    :param fastqs: list of (remote or local) fastq paths
    :type fastqs: list
    :param temp_dir: temporary directory
    :type temp_dir: str

    :return: all remote paths substituted with a local path
    :rtype: list
    """
    return [
        stream_file(fastq, os.path.join(temp_dir, os.path.basename(fastq)))
        if urlparse(fastq).scheme in ('http', 'https', 'ftp', 'ftps') else fastq
        for fastq in fastqs
    ]


def copy_or_create_whitelist(technology, bus_path, out_dir):
    """Copies a pre-packaged whitelist if it is provided. Otherwise, runs
    `bustools whitelist` to generate a whitelist.

    :param technology: single-cell technology used
    :type technology: str
    :param bus_path: path to BUS file generate the whitelist from
    :type bus_path: str
    :param out_dir: path to output directory
    :type out_dir: str

    :return: path to copied or generated whitelist
    :rtype: str
    """
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


def convert_transcripts_to_genes(txnames_path, t2g_path, genes_path):
    """Convert a textfile containing transcript IDs to another textfile containing
    gene IDs, given a transcript-to-gene mapping.

    :param txnames_path: path to transcripts.txt
    :type txnames_path: str
    :param t2g_path: path to transcript-to-genes mapping
    :type t2g_path: str
    :param genes_path: path to output genes.txt
    :type genes_path: str

    :return: path to written genes.txt
    :rtype: str
    """
    t2g = read_t2g(t2g_path)
    with open_as_text(txnames_path, 'r') as f, open_as_text(genes_path,
                                                            'w') as out:
        for line in f:
            if line.isspace():
                continue
            transcript = line.strip()
            if transcript not in t2g:
                logger.warning(
                    f'Transcript {transcript} was found in {txnames_path} but not in {t2g_path}. '
                    'This transcript will not be converted to a gene.'
                )
            attributes = t2g.get(transcript)

            if attributes:
                out.write(f'{attributes[0]}\n')
            else:
                out.write(f'{transcript}\n')
    return genes_path


@logger.namespaced('count')
def count(
    index_paths,
    t2g_path,
    technology,
    out_dir,
    fastqs,
    whitelist_path=None,
    tcc=False,
    mm=False,
    filter=None,
    kite=False,
    FB=False,
    temp_dir='tmp',
    threads=8,
    memory='4G',
    overwrite=False,
    loom=False,
    h5ad=False,
    cellranger=False,
    inspect=True,
    report=False,
):
    """Generates count matrices for single-cell RNA seq.

    :param index_paths: paths to kallisto indices
    :type index_paths: list
    :param t2g_path: path to transcript-to-gene mapping
    :type t2g_path: str
    :param technology: single-cell technology used
    :type technology: str
    :param out_dir: path to output directory
    :type out_dir: str
    :param fastqs: list of FASTQ file paths
    :type fastqs: list
    :param whitelist_path: path to whitelist, defaults to `None`
    :type whitelist_path: str, optional
    :param tcc: whether to generate a TCC matrix instead of a gene count matrix,
                defaults to `False`
    :type tcc: bool, optional
    :param mm: whether to include BUS records that pseudoalign to multiple genes,
               defaults to `False`
    :type mm: bool, optional
    :param filter: filter to use to generate a filtered count matrix,
                   defaults to `None`
    :type filter: str, optional
    :param kite: Whether this is a KITE workflow
    :type kite: bool, optional
    :param FB: whether 10x Genomics Feature Barcoding technology was used,
               defaults to `False`
    :type FB: bool, optional
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional
    :param threads: number of threads to use, defaults to `8`
    :type threads: int, optional
    :param memory: amount of memory to use, defaults to `4G`
    :type memory: str, optional
    :param overwrite: overwrite an existing index file, defaults to `False`
    :type overwrite: bool, optional
    :param loom: whether to convert the final count matrix into a loom file,
                 defaults to `False`
    :type loom: bool, optional
    :param h5ad: whether to convert the final count matrix into a h5ad file,
                 defaults to `False`
    :type h5ad: bool, optional
    :param cellranger: whether to convert the final count matrix into a
                       cellranger-compatible matrix, defaults to `False`
    :type cellranger: bool, optional
    :param inspect: whether or not to inspect the output BUS file and generate
                    the inspect.json
    :type inspect: bool, optional
    :param report: generate an HTMl report, defaults to `False`
    :type report: bool, optional

    :return: dictionary containing path to generated index
    :rtype: dict
    """
    STATS.start()
    if not isinstance(index_paths, list):
        index_paths = [index_paths]

    results = {}

    make_directory(out_dir)
    unfiltered_results = results.setdefault('unfiltered', {})

    bus_result = {
        'bus': os.path.join(out_dir, BUS_FILENAME),
        'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
        'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
        'info': os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    }
    if any(not os.path.exists(path)
           for name, path in bus_result.items()) or overwrite:
        if len(index_paths) > 1:
            bus_result = kallisto_bus_split(
                fastqs,
                index_paths,
                technology,
                out_dir,
                temp_dir=temp_dir,
                threads=threads,
                memory=memory
            )
        else:
            # Pipe any remote files.
            fastqs = stream_fastqs(fastqs, temp_dir=temp_dir)
            bus_result = kallisto_bus(
                fastqs, index_paths[0], technology, out_dir, threads=threads
            )
    else:
        logger.info(
            'Skipping kallisto bus because output files already exist. Use the --overwrite flag to overwrite.'
        )
    unfiltered_results.update(bus_result)

    sort_result = bustools_sort(
        bus_result['bus'],
        os.path.join(
            temp_dir,
            update_filename(os.path.basename(bus_result['bus']), SORT_CODE)
        ),
        temp_dir=temp_dir,
        threads=threads,
        memory=memory
    )
    if not whitelist_path:
        logger.info('Whitelist not provided')
        whitelist_path = copy_or_create_whitelist(
            technology, sort_result['bus'], out_dir
        )
        unfiltered_results.update({'whitelist': whitelist_path})

    prev_result = sort_result
    if FB:
        logger.info(f'Copying {technology} feature-to-barcode map to {out_dir}')
        map_path = copy_map(technology, out_dir)
        project_result = bustools_project(
            sort_result['bus'],
            os.path.join(
                temp_dir,
                update_filename(
                    os.path.basename(sort_result['bus']), PROJECT_CODE
                )
            ), map_path, bus_result['ecmap'], bus_result['txnames']
        )

        sort2_result = bustools_sort(
            project_result['bus'],
            os.path.join(
                temp_dir,
                update_filename(
                    os.path.basename(project_result['bus']), SORT_CODE
                )
            ),
            temp_dir=temp_dir,
            threads=threads,
            memory=memory
        )
        prev_result = sort2_result

    if inspect:
        inspect_result = bustools_inspect(
            prev_result['bus'], os.path.join(out_dir, INSPECT_FILENAME),
            whitelist_path, bus_result['ecmap']
        )
        unfiltered_results.update(inspect_result)
    correct_result = bustools_correct(
        prev_result['bus'],
        os.path.join(
            temp_dir,
            update_filename(os.path.basename(prev_result['bus']), CORRECT_CODE)
        ), whitelist_path
    )
    sort3_result = bustools_sort(
        correct_result['bus'],
        os.path.join(out_dir, f'output.{UNFILTERED_CODE}.bus'),
        temp_dir=temp_dir,
        threads=threads,
        memory=memory
    )
    unfiltered_results.update({'bus_scs': sort3_result['bus']})

    counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    make_directory(counts_dir)
    counts_prefix = os.path.join(
        counts_dir,
        TCC_PREFIX if tcc else FEATURE_PREFIX if kite else COUNTS_PREFIX
    )
    count_result = bustools_count(
        sort3_result['bus'],
        counts_prefix,
        t2g_path,
        bus_result['ecmap'],
        bus_result['txnames'],
        tcc=tcc,
        mm=mm,
    )
    unfiltered_results.update(count_result)

    # Convert outputs.
    if loom or h5ad:
        unfiltered_results.update(
            convert_matrix(
                counts_dir,
                count_result['mtx'],
                count_result['barcodes'],
                genes_path=count_result.get('genes'),
                t2g_path=t2g_path,
                ec_path=count_result.get('ec'),
                txnames_path=bus_result['txnames'],
                name=FEATURE_NAME if kite else GENE_NAME,
                loom=loom,
                h5ad=h5ad,
                tcc=tcc,
                threads=threads
            )
        )
    if cellranger:
        if not tcc:
            cr_result = matrix_to_cellranger(
                count_result['mtx'], count_result['barcodes'],
                count_result['genes'], t2g_path,
                os.path.join(counts_dir, CELLRANGER_DIR)
            )
            unfiltered_results.update({'cellranger': cr_result})
        else:
            logger.warning(
                'TCC matrices can not be converted to cellranger-compatible format.'
            )

    if filter == 'bustools':
        filtered_counts_prefix = os.path.join(
            out_dir, FILTERED_COUNTS_DIR,
            TCC_PREFIX if tcc else FEATURE_PREFIX if kite else COUNTS_PREFIX
        )
        filtered_whitelist_path = os.path.join(
            out_dir, FILTER_WHITELIST_FILENAME
        )
        filtered_bus_path = os.path.join(out_dir, f'output.{FILTERED_CODE}.bus')
        results['filtered'] = filter_with_bustools(
            sort3_result['bus'],
            bus_result['ecmap'],
            bus_result['txnames'],
            t2g_path,
            filtered_whitelist_path,
            filtered_bus_path,
            counts_prefix=filtered_counts_prefix,
            kite=kite,
            tcc=tcc,
            temp_dir=temp_dir,
            threads=threads,
            memory=memory,
            loom=loom,
            h5ad=h5ad
        )

    # Generate report.
    STATS.end()
    stats_path = STATS.save(os.path.join(out_dir, KB_INFO_FILENAME))
    results.update({'stats': stats_path})
    if report:
        nb_path = os.path.join(out_dir, REPORT_NOTEBOOK_FILENAME)
        html_path = os.path.join(out_dir, REPORT_HTML_FILENAME)
        logger.info(
            f'Writing report Jupyter notebook at {nb_path} and rendering it to {html_path}'
        )
        report_result = render_report(
            stats_path,
            bus_result['info'],
            inspect_result['inspect'],
            nb_path,
            html_path,
            count_result['mtx'],
            count_result.get('barcodes'),
            count_result.get('genes'),
            t2g_path,
            temp_dir=temp_dir
        )
        unfiltered_results.update(report_result)

        if tcc:
            logger.warning(
                'Plots for TCC matrices have not yet been implemented. The HTML report will not contain any plots.'
            )

    return results


@logger.namespaced('count_smartseq')
def count_smartseq(
    index_paths,
    t2g_path,
    technology,
    out_dir,
    fastq_pairs,
    cell_ids=None,
    temp_dir='tmp',
    threads=8,
    memory='4G',
    overwrite=False,
    loom=False,
    h5ad=False,
):
    """Generates gene or isoform count matrices from Smart-seq reads.
    """
    STATS.start()
    if not isinstance(index_paths, list):
        index_paths = [index_paths]

    # Smart-seq does not support multiple indices.
    if len(index_paths) > 1:
        raise Exception(
            f'Technology {technology} does not support multiple indices.'
        )

    # Smart-seq does not support fastq streaming.
    if any(urlparse(fastq).scheme in ('http', 'https', 'ftp', 'ftps')
           for pair in fastq_pairs
           for fastq in pair):
        raise Exception(
            f'Technology {technology} does not support FASTQ streaming.'
        )

    results = {}

    make_directory(out_dir)

    pseudo_result = {
        'mtx': os.path.join(out_dir, ABUNDANCE_FILENAME),
        'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
        'cells': os.path.join(out_dir, CELLS_FILENAME),
        'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
        'info': os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    }

    if any(not os.path.exists(path)
           for name, path in pseudo_result.items()) or overwrite:
        # Write batch information.
        batch_result = write_smartseq_batch(
            fastq_pairs,
            cell_ids if cell_ids else list(range(len(fastq_pairs))),
            os.path.join(out_dir, BATCH_FILENAME),
        )
        results.update(batch_result)

        pseudo_result = kallisto_pseudo(
            batch_result['batch'], index_paths[0], out_dir, threads=threads
        )
    else:
        logger.info(
            'Skipping kallisto pseudo because output files already exist. Use the --overwrite flag to overwrite.'
        )
    results.update(pseudo_result)

    # Manually write genes.txt
    # NOTE: there will be duplicated genes
    genes_path = os.path.join(out_dir, GENES_FILENAME)
    results['genes'] = convert_transcripts_to_genes(
        pseudo_result['txnames'], t2g_path, genes_path
    )

    # Convert outputs.
    if loom or h5ad:
        results.update(
            convert_matrix(
                out_dir,
                pseudo_result['mtx'],
                pseudo_result['cells'],
                results['genes'],
                t2g_path=t2g_path,
                loom=loom,
                h5ad=h5ad,
                threads=threads
            )
        )

    # Generate report.
    STATS.end()
    stats_path = STATS.save(os.path.join(out_dir, KB_INFO_FILENAME))
    results.update({'stats': stats_path})

    return results


@logger.namespaced('count_lamanno')
def count_velocity(
    index_paths,
    t2g_path,
    cdna_t2c_path,
    intron_t2c_path,
    technology,
    out_dir,
    fastqs,
    whitelist_path=None,
    tcc=False,
    mm=False,
    filter=None,
    temp_dir='tmp',
    threads=8,
    memory='4G',
    overwrite=False,
    loom=False,
    h5ad=False,
    cellranger=False,
    report=False,
    inspect=True,
    nucleus=False,
):
    """Generates RNA velocity matrices for single-cell RNA seq.

    :param index_paths: paths to kallisto indices
    :type index_paths: list
    :param t2g_path: path to transcript-to-gene mapping
    :type t2g_path: str
    :param cdna_t2c_path: path to cDNA transcripts-to-capture file
    :type cdna_t2c_path: str
    :param intron_t2c_path: path to intron transcripts-to-capture file
    :type intron_t2c_path: str
    :param technology: single-cell technology used
    :type technology: str
    :param out_dir: path to output directory
    :type out_dir: str
    :param fastqs: list of FASTQ file paths
    :type fastqs: list
    :param whitelist_path: path to whitelist, defaults to `None`
    :type whitelist_path: str, optional
    :param tcc: whether to generate a TCC matrix instead of a gene count matrix,
                defaults to `False`
    :type tcc: bool, optional
    :param mm: whether to include BUS records that pseudoalign to multiple genes,
               defaults to `False`
    :type mm: bool, optional
    :param filter: filter to use to generate a filtered count matrix,
                   defaults to `None`
    :type filter: str, optional
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional
    :param threads: number of threads to use, defaults to `8`
    :type threads: int, optional
    :param memory: amount of memory to use, defaults to `4G`
    :type memory: str, optional
    :param overwrite: overwrite an existing index file, defaults to `False`
    :type overwrite: bool, optional
    :param loom: whether to convert the final count matrix into a loom file,
                 defaults to `False`
    :type loom: bool, optional
    :param h5ad: whether to convert the final count matrix into a h5ad file,
                 defaults to `False`
    :type h5ad: bool, optional
    :param cellranger: whether to convert the final count matrix into a
                       cellranger-compatible matrix, defaults to `False`
    :type cellranger: bool, optional
    :param report: generate HTML reports, defaults to `False`
    :type report: bool, optional
    :param inspect: whether or not to inspect the output BUS file and generate
                    the inspect.json
    :type inspect: bool, optional
    :param nucleus: whether this is a single-nucleus experiment. if `True`, the
                    spliced and unspliced count matrices will be summed,
                    defaults to `False`
    :type nucleus: bool, optional

    :return: dictionary containing path to generated index
    :rtype: dict
    """
    STATS.start()
    if not isinstance(index_paths, list):
        index_paths = [index_paths]

    results = {}
    make_directory(out_dir)
    unfiltered_results = results.setdefault('unfiltered', {})

    bus_result = {
        'bus': os.path.join(out_dir, BUS_FILENAME),
        'ecmap': os.path.join(out_dir, ECMAP_FILENAME),
        'txnames': os.path.join(out_dir, TXNAMES_FILENAME),
        'info': os.path.join(out_dir, KALLISTO_INFO_FILENAME)
    }
    if any(not os.path.exists(path)
           for name, path in bus_result.items()) or overwrite:
        if len(index_paths) > 1:
            bus_result = kallisto_bus_split(
                fastqs,
                index_paths,
                technology,
                out_dir,
                temp_dir=temp_dir,
                threads=threads,
                memory=memory
            )
        else:
            # Pipe any remote files.
            fastqs = stream_fastqs(fastqs, temp_dir=temp_dir)
            bus_result = kallisto_bus(
                fastqs, index_paths[0], technology, out_dir, threads=threads
            )
    else:
        logger.info(
            'Skipping kallisto bus because output files already exist. Use the --overwrite flag to overwrite.'
        )
    unfiltered_results.update(bus_result)

    sort_result = bustools_sort(
        bus_result['bus'],
        os.path.join(
            temp_dir,
            update_filename(os.path.basename(bus_result['bus']), SORT_CODE)
        ),
        temp_dir=temp_dir,
        threads=threads,
        memory=memory
    )
    if not whitelist_path:
        logger.info('Whitelist not provided')
        whitelist_path = copy_or_create_whitelist(
            technology, sort_result['bus'], out_dir
        )
        unfiltered_results.update({'whitelist': whitelist_path})

    if inspect:
        inspect_result = bustools_inspect(
            sort_result['bus'], os.path.join(out_dir, INSPECT_FILENAME),
            whitelist_path, bus_result['ecmap']
        )
        unfiltered_results.update(inspect_result)
    correct_result = bustools_correct(
        sort_result['bus'],
        os.path.join(
            temp_dir,
            update_filename(os.path.basename(sort_result['bus']), CORRECT_CODE)
        ), whitelist_path
    )
    sort2_result = bustools_sort(
        correct_result['bus'],
        os.path.join(out_dir, f'output.{UNFILTERED_CODE}.bus'),
        temp_dir=temp_dir,
        threads=threads,
        memory=memory
    )
    unfiltered_results.update({'bus_scs': sort2_result['bus']})

    prefixes = [BUS_CDNA_PREFIX, BUS_INTRON_PREFIX]
    # The prefix and t2cs are swapped because we call bustools capture with
    # the --complement flag.
    prefix_to_t2c = {
        BUS_CDNA_PREFIX: intron_t2c_path,
        BUS_INTRON_PREFIX: cdna_t2c_path,
    }
    counts_dir = os.path.join(out_dir, UNFILTERED_COUNTS_DIR)
    make_directory(counts_dir)
    for prefix, t2c_path in prefix_to_t2c.items():
        capture_result = bustools_capture(
            sort2_result['bus'],
            os.path.join(temp_dir, '{}.bus'.format(prefix)), t2c_path,
            bus_result['ecmap'], bus_result['txnames']
        )
        sort_result = bustools_sort(
            capture_result['bus'],
            os.path.join(out_dir, f'{prefix}.{UNFILTERED_CODE}.bus'),
            temp_dir=temp_dir,
            threads=threads,
            memory=memory
        )

        if prefix not in unfiltered_results:
            unfiltered_results[prefix] = {}
        unfiltered_results[prefix].update(sort_result)

        if inspect:
            inspect_result = bustools_inspect(
                sort_result['bus'],
                os.path.join(
                    out_dir, update_filename(INSPECT_FILENAME, prefix)
                ), whitelist_path, bus_result['ecmap']
            )
            unfiltered_results[prefix].update(inspect_result)

        counts_prefix = os.path.join(counts_dir, prefix)
        count_result = bustools_count(
            sort_result['bus'],
            counts_prefix,
            t2g_path,
            bus_result['ecmap'],
            bus_result['txnames'],
            tcc=tcc,
            mm=mm,
        )
        unfiltered_results[prefix].update(count_result)

        if cellranger:
            if not tcc:
                cr_result = matrix_to_cellranger(
                    count_result['mtx'], count_result['barcodes'],
                    count_result['genes'], t2g_path,
                    os.path.join(counts_dir, f'{CELLRANGER_DIR}_{prefix}')
                )
                unfiltered_results[prefix].update({'cellranger': cr_result})
            else:
                logger.warning(
                    'TCC matrices can not be converted to cellranger-compatible format.'
                )

    if loom or h5ad:
        unfiltered_results.update(
            convert_matrices(
                counts_dir,
                [unfiltered_results[prefix]['mtx'] for prefix in prefixes],
                [unfiltered_results[prefix]['barcodes'] for prefix in prefixes],
                genes_paths=[
                    unfiltered_results[prefix].get('genes')
                    for prefix in prefixes
                ],
                t2g_path=t2g_path,
                ec_paths=[
                    unfiltered_results[prefix].get('ec') for prefix in prefixes
                ],
                txnames_path=bus_result['txnames'],
                loom=loom,
                h5ad=h5ad,
                tcc=tcc,
                nucleus=nucleus
            )
        )

    if filter:
        filtered_results = results.setdefault('filtered', {})
        if filter == 'bustools':
            filtered_results.update(
                filter_with_bustools(
                    sort2_result['bus'],
                    bus_result['ecmap'],
                    bus_result['txnames'],
                    t2g_path,
                    os.path.join(out_dir, FILTER_WHITELIST_FILENAME),
                    os.path.join(out_dir, f'output.{FILTERED_CODE}.bus'),
                    temp_dir=temp_dir,
                    memory=memory,
                    count=False
                )
            )

            for prefix, t2c_path in prefix_to_t2c.items():
                filtered_capture_result = bustools_capture(
                    filtered_results['bus_scs'],
                    os.path.join(temp_dir, '{}.bus'.format(prefix)), t2c_path,
                    bus_result['ecmap'], bus_result['txnames']
                )
                filtered_sort_result = bustools_sort(
                    filtered_capture_result['bus'],
                    os.path.join(out_dir, f'{prefix}.{FILTERED_CODE}.bus'),
                    temp_dir=temp_dir,
                    threads=threads,
                    memory=memory
                )

                filtered_results.setdefault(prefix,
                                            {}).update(filtered_sort_result)

                filtered_counts_dir = os.path.join(out_dir, FILTERED_COUNTS_DIR)
                make_directory(filtered_counts_dir)
                filtered_counts_prefix = os.path.join(
                    filtered_counts_dir, prefix
                )
                count_result = bustools_count(
                    filtered_sort_result['bus'],
                    filtered_counts_prefix,
                    t2g_path,
                    bus_result['ecmap'],
                    bus_result['txnames'],
                    tcc=tcc,
                    mm=mm,
                )
                filtered_results[prefix].update(count_result)

                if cellranger:
                    if not tcc:
                        cr_result = matrix_to_cellranger(
                            count_result['mtx'], count_result['barcodes'],
                            count_result['genes'], t2g_path,
                            os.path.join(
                                filtered_counts_dir,
                                f'{CELLRANGER_DIR}_{prefix}'
                            )
                        )
                        unfiltered_results[prefix].update({
                            'cellranger': cr_result
                        })
                    else:
                        logger.warning(
                            'TCC matrices can not be converted to cellranger-compatible format.'
                        )

        if loom or h5ad:
            filtered_results.update(
                convert_matrices(
                    filtered_counts_dir,
                    [filtered_results[prefix]['mtx'] for prefix in prefixes],
                    [
                        filtered_results[prefix]['barcodes']
                        for prefix in prefixes
                    ],
                    genes_paths=[
                        filtered_results[prefix].get('genes')
                        for prefix in prefixes
                    ],
                    t2g_path=t2g_path,
                    ec_paths=[
                        filtered_results[prefix].get('ec')
                        for prefix in prefixes
                    ],
                    txnames_path=bus_result['txnames'],
                    loom=loom,
                    h5ad=h5ad,
                    tcc=tcc,
                    nucleus=nucleus,
                )
            )

    STATS.end()
    stats_path = STATS.save(os.path.join(out_dir, KB_INFO_FILENAME))
    results.update({'stats': stats_path})

    # Reports
    nb_path = os.path.join(out_dir, REPORT_NOTEBOOK_FILENAME)
    html_path = os.path.join(out_dir, REPORT_HTML_FILENAME)
    if report:
        logger.info(
            f'Writing report Jupyter notebook at {nb_path} and rendering it to {html_path}'
        )
        report_result = render_report(
            stats_path,
            bus_result['info'],
            unfiltered_results['inspect'],
            nb_path,
            html_path,
            temp_dir=temp_dir
        )

        unfiltered_results.update(report_result)

        for prefix in prefix_to_t2c:
            nb_path = os.path.join(
                out_dir, update_filename(REPORT_NOTEBOOK_FILENAME, prefix)
            )
            html_path = os.path.join(
                out_dir, update_filename(REPORT_HTML_FILENAME, prefix)
            )
            logger.info(
                f'Writing report Jupyter notebook at {nb_path} and rendering it to {html_path}'
            )
            report_result = render_report(
                stats_path,
                bus_result['info'],
                unfiltered_results[prefix]['inspect'],
                nb_path,
                html_path,
                unfiltered_results[prefix]['mtx'],
                unfiltered_results[prefix].get('barcodes'),
                unfiltered_results[prefix].get('genes'),
                t2g_path,
                temp_dir=temp_dir
            )
            unfiltered_results[prefix].update(report_result)
        if tcc:
            logger.warning(
                'Plots for TCC matrices have not yet been implemented. The HTML report will not contain any plots.'
            )

    return results

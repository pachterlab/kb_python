import glob
import logging
import os
import tarfile

from .config import get_kallisto_binary_path
from .fasta import (
    FASTA,
    generate_cdna_fasta,
    generate_intron_fasta,
    generate_kite_fasta,
)
from .gtf import GTF
from .utils import (
    concatenate_files,
    decompress_gzip,
    download_file,
    get_temporary_filename,
    open_as_text,
    run_executable,
)

logger = logging.getLogger(__name__)


def sort_gtf(gtf_path, out_path):
    """Sorts a GTF file based on its chromosome, start position, line number.

    :param gtf_path: path to GTF file
    :type gtf_path: str

    :return: path to sorted GTF file, set of chromosomes in GTF file
    :rtype: tuple
    """
    logger.info(f'Sorting {gtf_path} to {out_path}')
    gtf = GTF(gtf_path)
    return gtf.sort(out_path)


def sort_fasta(fasta_path, out_path):
    """Sorts a FASTA file based on its header.

    :param fasta_path: path to FASTA file
    :type fasta_path: str

    :return: path to sorted FASTA file, set of chromosomes in FASTA file
    :rtype: tuple
    """
    logger.info(f'Sorting {fasta_path} to {out_path}')
    fasta = FASTA(fasta_path)
    return fasta.sort(out_path)


def check_chromosomes(fasta_chromosomes, gtf_chromosomes):
    """Compares the two chromosome sets and outputs warnings if there are
    unique chromosomes in either set.

    :param fasta_chromosomes: set of chromosomes found in FASTA
    :type fasta_chromosomes: set
    :param gtf_chromosomes: set of chromosomes found in GTF
    :type gtf_chromosomes: set

    :return: intersection of the two sets
    :rtype: set
    """
    fasta_unique = fasta_chromosomes - gtf_chromosomes
    gtf_unique = gtf_chromosomes - fasta_chromosomes
    if fasta_unique:
        logger.warning((
            'The following chromosomes were found in the FASTA but does not have '
            'any "transcript" features in the GTF: {}. '
            'No sequences will be generated for these chromosomes.'
        ).format(', '.join(fasta_unique)))
    if gtf_unique:
        logger.warning((
            'The following chromosomes were found to have "transcript" features '
            'in the GTF but does not exist in the FASTA. '
            'No sequences will be generated for these chromosomes.'
        ))
    chromosomes = set.intersection(fasta_chromosomes, gtf_chromosomes)

    return chromosomes


def create_t2g_from_fasta(fasta_path, t2g_path):
    """Parse FASTA headers to get transcripts-to-gene mapping.

    :param fasta_path: path to FASTA file
    :type fasta_path: str
    :param t2g_path: path to output transcript-to-gene mapping
    :type t2g_path: str

    :return: dictionary containing path to generated t2g mapping
    :rtype: dict
    """
    logger.info('Creating transcript-to-gene mapping at {}'.format(t2g_path))
    with open_as_text(t2g_path, 'w') as f:
        fasta = FASTA(fasta_path)
        for info, _ in fasta.entries():
            if 'feature_id' in info['group']:
                row = [
                    info['sequence_id'],
                    info['group']['feature_id'],
                    info['group']['feature_id'],
                ]
            else:
                row = [
                    info['sequence_id'], info['group']['gene_id'],
                    info['group'].get('gene_name', ''),
                    info['group'].get('transcript_name',
                                      ''), info['group'].get('chr', ''),
                    info['group'].get('start',
                                      ''), info['group'].get('end', ''),
                    info['group'].get('strand', '')
                ]
            f.write('\t'.join(str(item) for item in row) + '\n')
    return {'t2g': t2g_path}


def create_t2g_from_gtf(gtf_path, t2g_path, intron=False):
    """Creates a transcript-to-gene mapping from a GTF file.

    GTF entries that have `transcript` as its feature are parsed for
    the `transcript_id`, `gene_id` and `gene_name`.

    :param gtf_path: path to GTF file
    :type gtf_path: str
    :param t2g_path: path to output transcript-to-gene mapping
    :type t2g_path: str
    :param intron: whether or not to include intron transcript ids (with the
                   `-I` prefix), defaults to `False`
    :type intron: bool, optional

    :return: dictionary containing path to generated t2g mapping
    :rtype: dict
    """
    logger.info('Creating transcript-to-gene mapping at {}'.format(t2g_path))
    gtf = GTF(gtf_path)
    with open_as_text(t2g_path, 'w') as f:
        for entry in gtf.entries():
            if entry['feature'] == 'transcript':
                transcript_id = entry['group']['transcript_id']
                transcript_version = entry['group'].get(
                    'transcript_version', None
                )
                transcript = '{}.{}'.format(
                    transcript_id, transcript_version
                ) if transcript_version else transcript_id
                gene_id = entry['group']['gene_id']
                gene_version = entry['group'].get('gene_version', None)
                gene = '{}.{}'.format(
                    gene_id, gene_version
                ) if gene_version else gene_id

                row = [
                    transcript,
                    gene,
                    entry['group'].get('gene_name', ''),
                    entry['group'].get('transcript_name', ''),
                    entry.get('seqname', ''),
                    entry.get('start', ''),
                    entry.get('end', ''),
                    entry.get('strand', ''),
                ]
                f.write('\t'.join(str(item) for item in row) + '\n')

                if intron:
                    intron_row = row.copy()
                    intron_row[0] = intron_row[0] + '-I'
                    f.write('\t'.join(str(item) for item in intron_row) + '\n')

    return {'t2g': t2g_path}


def create_t2c(fasta_path, t2c_path):
    """Creates a transcripts-to-capture list from a FASTA file.

    :param fasta_path: path to FASTA file
    :type fasta_path: str
    :param t2c_path: path to output transcripts-to-capture list
    :type t2c_path: str

    :return: dictionary containing path to generated t2c list
    :rtype: dict
    """
    fasta = FASTA(fasta_path)
    with open_as_text(t2c_path, 'w') as f:
        for info, _ in fasta.entries():
            sequence_id = info['sequence_id']
            f.write('{}\n'.format(sequence_id))
    return {'t2c': t2c_path}


def kallisto_index(fasta_path, index_path, k=31):
    """Runs `kallisto index`.

    :param fasta_path: path to FASTA file
    :type fasta_path: str
    :param index_path: path to output kallisto index
    :type index_path: str
    :param k: k-mer length, defaults to 31
    :type k: int, optional

    :return: dictionary containing path to generated index
    :rtype: dict
    """
    logger.info(f'Indexing {fasta_path} to {index_path}')
    command = [
        get_kallisto_binary_path(), 'index', '-i', index_path, '-k', k,
        fasta_path
    ]
    run_executable(command)
    return {'index': index_path}


def split_and_index(fasta_path, index_prefix, n=2, k=31, temp_dir='tmp'):
    """Split a FASTA file into `n` parts and index each one.

    :param fasta_path: path to FASTA file
    :type fasta_path: str
    :param index_prefix: prefix of output kallisto indices
    :type index_prefix: str
    :param n: split the index into `n` files, defaults to `2`
    :type n: int, optional
    :param k: k-mer length, defaults to 31
    :type k: int, optional
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional

    :return: dictionary containing path to generated index
    :rtype: dict
    """
    fastas = []
    indices = []

    logger.info(f'Splitting {fasta_path} into {n} parts')
    size = int(os.path.getsize(fasta_path) / n) + 4

    fasta = FASTA(fasta_path)
    fasta_entries = fasta.entries(parse=False)
    finished = False
    for i in range(n):
        fasta_part_path = get_temporary_filename(temp_dir)
        index_part_path = f'{index_prefix}.{i}'
        fastas.append(fasta_part_path)
        indices.append(index_part_path)

        with open(fasta_part_path, 'w') as f:
            logger.debug(f'Writing {fasta_part_path}')
            while f.tell() < size:
                try:
                    info, seq = next(fasta_entries)
                except StopIteration:
                    finished = True
                    break
                f.write(f'{info}\n{seq}\n')

        if finished:
            break

    built = []
    for fasta_part_path, index_part_path in zip(fastas, indices):
        result = kallisto_index(fasta_part_path, index_part_path, k=k)
        built.append(result['index'])

    return {'indices': built}


def download_reference(reference, files, temp_dir='tmp', overwrite=False):
    """Downloads a provided reference file from a static url.

    The configuration for provided references is in `config.py`.

    :param reference: a Reference object, as defined in `config.py`
    :type reference: Reference
    :param files: dictionary that has the command-line option as keys and
                  the path as values. used to determine if all the required
                  paths to download the given reference have been provided
    :type files: dict
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional
    :param overwrite: overwrite an existing index file, defaults to `False`
    :type overwrite: bool, optional

    :raises Exception: if the required options are not provided

    :return: dictionary containing paths to generated file(s)
    :rtype: dict
    """
    results = {}
    if not any(os.path.exists(file) for file in files) or overwrite:
        # Make sure all the required file paths are there.
        diff = set(reference.files.keys()) - set(files.keys())
        if diff:
            raise Exception(
                'the following options are required to download this reference: {}'
                .format(','.join(diff))
            )

        url = reference.url
        path = os.path.join(temp_dir, os.path.basename(url))
        logging.info(
            'Downloading files for {} from {} to {}'.format(
                reference.name, url, path
            )
        )
        local_path = download_file(url, path)

        logging.info('Extracting files from {}'.format(local_path))
        with tarfile.open(local_path, 'r:gz') as f:
            f.extractall(temp_dir)

        for option in reference.files:
            os.rename(
                os.path.join(temp_dir, reference.files[option]), files[option]
            )
            results.update({option: files[option]})
    else:
        logger.info(
            'Skipping download because some files already exist. Use the --overwrite flag to overwrite.'
        )
    return results


def decompress_file(path, temp_dir='tmp'):
    """Decompress the given path if it is a .gz file. Otherwise, return the
    original path.

    :param path: path to the file
    :type path: str

    :return: unaltered `path` if the file is not a .gz file, otherwise path to the
             uncompressed file
    :rtype: str
    """
    if path.endswith('.gz'):
        logger.info('Decompressing {} to {}'.format(path, temp_dir))
        return decompress_gzip(
            path,
            os.path.join(temp_dir,
                         os.path.splitext(os.path.basename(path))[0])
        )
    else:
        return path


def ref(
    fasta_paths,
    gtf_paths,
    cdna_path,
    index_path,
    t2g_path,
    n=1,
    k=None,
    temp_dir='tmp',
    overwrite=False
):
    """Generates files necessary to generate count matrices for single-cell RNA-seq.

    :param fasta_paths: list of paths to genomic FASTA files
    :type fasta_paths: list
    :param gtf_paths: list of paths to GTF files
    :type gtf_paths: list
    :param cdna_path: path to generate the cDNA FASTA file
    :type cdna_path: str
    :param t2g_path: path to output transcript-to-gene mapping
    :type t2g_path: str
    :param n: split the index into `n` files
    :type n: int
    :param k: override default kmer length 31, defaults to `None`
    :type k: int, optional
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional
    :param overwrite: overwrite an existing index file, defaults to `False`
    :type overwrite: bool, optional

    :return: dictionary containing paths to generated file(s)
    :rtype: dict
    """
    if not isinstance(fasta_paths, list):
        fasta_paths = [fasta_paths]
    if not isinstance(gtf_paths, list):
        gtf_paths = [gtf_paths]

    results = {}
    t2gs = []
    cdnas = []
    for fasta_path, gtf_path in zip(fasta_paths, gtf_paths):
        logger.info(f'Preparing {fasta_path}, {gtf_path}')
        gtf_path = decompress_file(gtf_path, temp_dir=temp_dir)
        t2g_result = create_t2g_from_gtf(
            gtf_path, get_temporary_filename(temp_dir)
        )
        t2gs.append(t2g_result['t2g'])

        if not glob.glob(f'{index_path}*') or overwrite:
            fasta_path = decompress_file(fasta_path, temp_dir=temp_dir)
            sorted_fasta_path, fasta_chromosomes = sort_fasta(
                fasta_path, get_temporary_filename(temp_dir)
            )
            sorted_gtf_path, gtf_chromosomes = sort_gtf(
                gtf_path, get_temporary_filename(temp_dir)
            )

            cdna_temp_path = get_temporary_filename(temp_dir)
            logger.info(
                f'Splitting genome {fasta_path} into cDNA at {cdna_temp_path}'
            )
            chromosomes = check_chromosomes(fasta_chromosomes, gtf_chromosomes)
            cdna_fasta_path = generate_cdna_fasta(
                sorted_fasta_path,
                sorted_gtf_path,
                cdna_temp_path,
                chromosomes=chromosomes
            )
            cdnas.append(cdna_fasta_path)

    # Concatenate t2gs
    logger.info(
        f'Concatenating {len(t2gs)} transcript-to-gene mappings to {t2g_path}'
    )
    t2g_path = concatenate_files(*t2gs, out_path=t2g_path, temp_dir=temp_dir)
    results.update({'t2g': t2g_path})

    # Concatenate cdnas and index
    if not glob.glob(f'{index_path}*') or overwrite:
        logger.info(f'Concatenating {len(cdnas)} cDNAs to {cdna_path}')
        cdna_fasta_path = concatenate_files(
            *cdnas, out_path=cdna_path, temp_dir=temp_dir
        )

        if k and k != 31:
            logger.warning(
                f'Using provided k-mer length {k} instead of optimal length 31'
            )
        index_result = split_and_index(
            cdna_fasta_path, index_path, n=n, k=k or 31, temp_dir=temp_dir
        ) if n > 1 else kallisto_index(
            cdna_fasta_path, index_path, k=k or 31
        )
        results.update(index_result)
    else:
        logger.info(
            'Skipping kallisto index because {} already exists. Use the --overwrite flag to overwrite.'
            .format(index_path)
        )

    return results


def ref_kite(
    feature_path,
    fasta_path,
    index_path,
    t2g_path,
    n=1,
    k=None,
    no_mismatches=False,
    temp_dir='tmp',
    overwrite=False
):
    """Generates files necessary for feature barcoding with the KITE workflow.

    :param feature_path: path to TSV containing barcodes and feature names
    :type feature_path: str
    :param fasta_path: path to generate fasta file containing all sequences
                       that are 1 hamming distance from the provide barcodes
                       (including the actual sequence)
    :type fasta_path: str
    :param t2g_path: path to output transcript-to-gene mapping
    :type t2g_path: str
    :param n: split the index into `n` files
    :type n: int
    :param k: override calculated optimal kmer length, defaults to `None`
    :type k: int, optional
    :param no_mismatches: whether to generate hamming distance 1 variants,
                          defaults to `False`
    :type no_mismatches: bool, optional
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional
    :param overwrite: overwrite an existing index file, defaults to `False`
    :type overwrite: bool, optional

    :return: dictionary containing paths to generated file(s)
    :rtype: dict
    """
    results = {}
    feature_path = decompress_file(feature_path, temp_dir=temp_dir)
    logger.info('Generating mismatch FASTA at {}'.format(fasta_path))
    kite_path, length = generate_kite_fasta(
        feature_path, fasta_path, no_mismatches=no_mismatches
    )
    results.update({'fasta': kite_path})
    t2g_result = create_t2g_from_fasta(fasta_path, t2g_path)
    results.update(t2g_result)

    if not glob.glob(f'{index_path}*') or overwrite:
        optimal_k = length if length % 2 else length - 1
        if k and k != optimal_k:
            logger.warning(
                f'Using provided k-mer length {k} instead of calculated optimal length {optimal_k}'
            )
        index_result = split_and_index(
            kite_path, index_path, n=n, k=k or optimal_k, temp_dir=temp_dir
        ) if n > 1 else kallisto_index(
            kite_path, index_path, k=k or optimal_k
        )
        results.update(index_result)
    else:
        logger.info(
            'Skipping kallisto index because {} already exists. Use the --overwrite flag to overwrite.'
            .format(index_path)
        )
    return results


def ref_lamanno(
    fasta_paths,
    gtf_paths,
    cdna_path,
    intron_path,
    index_path,
    t2g_path,
    cdna_t2c_path,
    intron_t2c_path,
    n=1,
    k=None,
    flank=None,
    temp_dir='tmp',
    overwrite=False,
):
    """Generates files necessary to generate RNA velocity matrices for single-cell RNA-seq.

    :param fasta_paths: list of paths to genomic FASTA files
    :type fasta_paths: list
    :param gtf_paths: list of paths to GTF files
    :type gtf_paths: list
    :param cdna_path: path to generate the cDNA FASTA file
    :type cdna_path: str
    :param intron_path: path to generate the intron FASTA file
    :type intron_path: str
    :param t2g_path: path to output transcript-to-gene mapping
    :type t2g_path: str
    :param cdna_t2c_path: path to generate the cDNA transcripts-to-capture file
    :type cdna_t2c_path: str
    :param intron_t2c_path: path to generate the intron transcripts-to-capture file
    :type intron_t2c_path: str
    :param n: split the index into `n` files
    :type n: int
    :param k: override default kmer length (31), defaults to `None`
    :type k: int, optional
    :param flank: number of bases to include from the flanking regions
                  when generating the intron FASTA, defaults to `None`, which
                  sets the flanking region to be k - 1 bases.
    :type flank: int, optional
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional
    :param overwrite: overwrite an existing index file, defaults to `False`
    :type overwrite: bool, optional

    :return: dictionary containing paths to generated file(s)
    :rtype: dict
    """
    if not isinstance(fasta_paths, list):
        fasta_paths = [fasta_paths]
    if not isinstance(gtf_paths, list):
        gtf_paths = [gtf_paths]

    results = {}
    cdnas = []
    introns = []
    cdna_t2cs = []
    intron_t2cs = []
    if not glob.glob(f'{index_path}*') or overwrite:
        for fasta_path, gtf_path in zip(fasta_paths, gtf_paths):
            logger.info(f'Preparing {fasta_path}, {gtf_path}')

            fasta_path = decompress_file(fasta_path, temp_dir=temp_dir)
            sorted_fasta_path, fasta_chromosomes = sort_fasta(
                fasta_path, get_temporary_filename(temp_dir)
            )
            gtf_path = decompress_file(gtf_path, temp_dir=temp_dir)
            sorted_gtf_path, gtf_chromosomes = sort_gtf(
                gtf_path, get_temporary_filename(temp_dir)
            )
            cdna_temp_path = get_temporary_filename(temp_dir)
            logger.info(
                f'Splitting genome {fasta_path} into cDNA at {cdna_temp_path}'
            )
            chromosomes = check_chromosomes(fasta_chromosomes, gtf_chromosomes)
            cdna_fasta_path = generate_cdna_fasta(
                sorted_fasta_path,
                sorted_gtf_path,
                cdna_temp_path,
                chromosomes=chromosomes
            )
            cdnas.append(cdna_fasta_path)
            cdna_t2c_temp_path = get_temporary_filename(temp_dir)
            logger.info(
                f'Creating cDNA transcripts-to-capture at {cdna_t2c_temp_path}'
            )
            cdna_t2c_result = create_t2c(cdna_fasta_path, cdna_t2c_temp_path)
            cdna_t2cs.append(cdna_t2c_result['t2c'])
            intron_temp_path = get_temporary_filename(temp_dir)
            logger.info(f'Splitting genome into introns at {intron_temp_path}')
            intron_fasta_path = generate_intron_fasta(
                sorted_fasta_path,
                sorted_gtf_path,
                intron_temp_path,
                chromosomes=chromosomes,
                flank=flank if flank is not None else k -
                1 if k is not None else 30
            )
            introns.append(intron_fasta_path)
            intron_t2c_temp_path = get_temporary_filename(temp_dir)
            logger.info(
                f'Creating intron transcripts-to-capture at {intron_t2c_temp_path}'
            )
            intron_t2c_result = create_t2c(
                intron_fasta_path, intron_t2c_temp_path
            )
            intron_t2cs.append(intron_t2c_result['t2c'])

        # Concatenate
        logger.info(f'Concatenating {len(cdnas)} cDNA FASTAs to {cdna_path}')
        cdna_path = concatenate_files(
            *cdnas, out_path=cdna_path, temp_dir=temp_dir
        )
        logger.info(
            f'Concatenating {len(cdna_t2cs)} cDNA transcripts-to-captures to {cdna_t2c_path}'
        )
        cdna_t2c_path = concatenate_files(
            *cdna_t2cs, out_path=cdna_t2c_path, temp_dir=temp_dir
        )
        logger.info(
            f'Concatenating {len(introns)} intron FASTAs to {intron_path}'
        )
        intron_path = concatenate_files(
            *introns, out_path=intron_path, temp_dir=temp_dir
        )
        logger.info(
            f'Concatenating {len(intron_t2cs)} intron transcripts-to-captures to {intron_t2c_path}'
        )
        intron_t2c_path = concatenate_files(
            *intron_t2cs, out_path=intron_t2c_path, temp_dir=temp_dir
        )
        results.update({
            'cdna_fasta': cdna_path,
            'cdna_t2c': cdna_t2c_path,
            'intron_fasta': intron_path,
            'intron_t2c': intron_t2c_path
        })

        # Concatenate cDNA and intron fastas to generate T2G and build index
        combined_path = get_temporary_filename(temp_dir)
        logger.info(f'Concatenating cDNA and intron FASTAs to {combined_path}')
        combined_path = concatenate_files(
            cdna_path, intron_path, out_path=combined_path, temp_dir=temp_dir
        )
        t2g_result = create_t2g_from_fasta(combined_path, t2g_path)
        results.update(t2g_result)
        if k and k != 31:
            logger.warning(
                f'Using provided k-mer length {k} instead of optimal length 31'
            )

        # If n = 1, make single index
        # if n = 2, make two indices, one for spliced and another for unspliced
        # if n > 2, make n indices, one for spliced, another n - 1 for unspliced
        if n == 1:
            index_result = kallisto_index(combined_path, index_path, k=k or 31)
        else:
            cdna_index_result = kallisto_index(
                cdna_path, f'{index_path}_cdna', k=k or 31
            )
            if n == 2:
                intron_index_result = kallisto_index(
                    intron_path, f'{index_path}_intron', k=k or 31
                )
                index_result = {
                    'indices': [
                        cdna_index_result['index'], intron_index_result['index']
                    ]
                }
            else:
                split_index_result = split_and_index(
                    intron_path,
                    f'{index_path}_intron',
                    n=n - 1,
                    k=k or 31,
                    temp_dir=temp_dir
                )
                index_result = {
                    'indices': [
                        cdna_index_result['index'],
                        *split_index_result['indices']
                    ]
                }
        results.update(index_result)
    else:
        logger.info(
            'Skipping kallisto index because {} already exists. Use the --overwrite flag to overwrite.'
            .format(index_path)
        )

    return results

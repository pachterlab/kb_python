import glob
import itertools
import os
import tarfile

import ngs_tools as ngs
import pandas as pd

from .config import get_kallisto_binary_path
from .logging import logger
from .utils import (
    concatenate_files,
    decompress_gzip,
    download_file,
    get_temporary_filename,
    open_as_text,
    run_executable,
)


def generate_kite_fasta(feature_path, out_path, no_mismatches=False):
    """Generate a FASTA file for feature barcoding with the KITE workflow.

    This FASTA contains all sequences that are 1 hamming distance from the
    provided barcodes. The file of barcodes must be a 2-column TSV containing
    the barcode sequences in the first column and their corresponding feature
    name in the second column. If hamming distance 1 variants collide for any
    pair of barcodes, the hamming distance 1 variants for those barcodes are
    not generated.

    :param feature_path: path to TSV containing barcodes and feature names
    :type feature_path: str
    :param out_path: path to FASTA to generate
    :type out_path: str
    :param no_mismatches: whether to generate hamming distance 1 variants,
                          defaults to `False`
    :type no_mismatches: bool, optional

    :raises Exception: if there are barcodes of different lengths
    :raises Exception: if there are duplicate barcodes

    :return: (path to generated FASTA, set of barcode lengths)
    :rtype: tuple
    """

    def generate_mismatches(name, sequence):
        """Helper function to generate 1 hamming distance mismatches.

        :param name: name of the sequence
        :type name: str
        :param sequence: the sequence
        :type sequence: str

        :return: a generator that yields a tuple of (unique name, mismatched sequence)
        :rtype: generator
        """
        sequence = sequence.upper()
        for i in range(len(sequence)):
            base = sequence[i]
            before = sequence[:i]
            after = sequence[i + 1:]

            for j, different in enumerate([b for b in ['A', 'C', 'G', 'T']
                                           if b != base]):
                yield f'{name}-{i}.{j+1}', f'{before}{different}{after}'

    df_features = pd.read_csv(
        feature_path, sep='\t', header=None, names=['sequence', 'name']
    )

    lengths = set()
    features = {}
    variants = {}
    # Generate all feature barcode variations before saving to check for collisions.
    for i, row in df_features.iterrows():
        # Check that the first column contains the sequence
        # and the second column the feature name.
        if ngs.sequence.SEQUENCE_PARSER.search(row.sequence.upper()):
            raise Exception((
                'Encountered non-ATCG basepairs in barcode sequence {row.sequence}. '
                'Does the first column contain the sequences and the second column the feature names?'
            ))

        lengths.add(len(row.sequence))
        features[row['name']] = row.sequence
        variants[row['name']] = {
            name: seq
            for name, seq in generate_mismatches(row['name'], row.sequence)
            if not no_mismatches
        }

    # Check duplicate barcodes.
    duplicates = set([
        bc for bc in features.values() if list(features.values()).count(bc) > 1
    ])
    if len(duplicates) > 0:
        raise Exception(
            'Duplicate feature barcodes: {}'.format(' '.join(duplicates))
        )
    if len(lengths) > 1:
        logger.warning(
            'Detected barcodes of different lengths: {}'.format(
                ','.join(str(l) for l in lengths)  # noqa
            )
        )
    for f1, f2 in itertools.combinations(variants.keys(), 2):
        v1 = variants[f1]
        v2 = variants[f2]
        if len(set.intersection(set(v1.values()), set(v2.values()))) > 0:
            logger.warning((
                f'Detected features with barcodes that are not >2 hamming distance apart: {f1}, {f2}. '
                f'Hamming distance 1 variants for these barcodes will not be generated.'
            ))
            variants[f1] = {}
            variants[f2] = {}

    # Write FASTA
    with ngs.fasta.Fasta(out_path, 'w') as f:
        for feature, barcode in features.items():
            attributes = {'feature_id': feature}
            header = ngs.fasta.FastaEntry.make_header(feature, attributes)
            entry = ngs.fasta.FastaEntry(header, barcode)
            f.write(entry)

            for name, variant in variants[feature].items():
                header = ngs.fasta.FastaEntry.make_header(name, attributes)
                entry = ngs.fasta.FastaEntry(header, variant)
                f.write(entry)

    return out_path, min(lengths)


def create_t2g_from_fasta(fasta_path, t2g_path):
    """Parse FASTA headers to get transcripts-to-gene mapping.

    :param fasta_path: path to FASTA file
    :type fasta_path: str
    :param t2g_path: path to output transcript-to-gene mapping
    :type t2g_path: str

    :return: dictionary containing path to generated t2g mapping
    :rtype: dict
    """
    logger.info(f'Creating transcript-to-gene mapping at {t2g_path}')
    with ngs.fasta.Fasta(fasta_path, 'r') as f_in, open_as_text(t2g_path,
                                                                'w') as f_out:
        for entry in f_in:
            attributes = entry.attributes

            if 'feature_id' in attributes:
                feature_id = attributes['feature_id']
                row = [entry.name, feature_id, feature_id]
            else:
                gene_id = attributes['gene_id']
                gene_name = attributes.get('gene_name', '')
                transcript_name = attributes.get('transcript_name', '')
                chromosome = attributes['chr']
                start = attributes['start']
                end = attributes['end']
                strand = attributes['strand']
                row = [
                    entry.name,
                    gene_id,
                    gene_name,
                    transcript_name,
                    chromosome,
                    start,
                    end,
                    strand,
                ]
            f_out.write('\t'.join(str(item) for item in row) + '\n')

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
    with ngs.fasta.Fasta(fasta_path, 'r') as f_in, open_as_text(t2c_path,
                                                                'w') as f_out:
        for entry in f_in:
            f_out.write(f'{entry.name}\n')
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

    with ngs.fasta.Fasta(fasta_path, 'r') as f_in:
        fasta_iter = iter(f_in)
        finished = False
        for i in range(n):
            fasta_part_path = get_temporary_filename(temp_dir)
            index_part_path = f'{index_prefix}.{i}'
            fastas.append(fasta_part_path)
            indices.append(index_part_path)

            with ngs.fasta.Fasta(fasta_part_path, 'w') as f_out:
                logger.debug(f'Writing {fasta_part_path}')
                while f_out.tell() < size:
                    try:
                        entry = next(fasta_iter)
                    except StopIteration:
                        finished = True
                        break
                    f_out.write(entry)

            if finished:
                break

    built = []
    for fasta_part_path, index_part_path in zip(fastas, indices):
        result = kallisto_index(fasta_part_path, index_part_path, k=k)
        built.append(result['index'])

    return {'indices': built}


@logger.namespaced('download')
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
    if not ngs.utils.all_exists(*list(files.values())) or overwrite:
        # Make sure all the required file paths are there.
        diff = set(reference.files.keys()) - set(files.keys())
        if diff:
            raise Exception(
                'the following options are required to download this reference: {}'
                .format(','.join(diff))
            )

        url = reference.url
        path = os.path.join(temp_dir, os.path.basename(url))
        logger.info(
            'Downloading files for {} from {} to {}'.format(
                reference.name, url, path
            )
        )
        local_path = download_file(url, path)

        logger.info('Extracting files from {}'.format(local_path))
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


@logger.namespaced('ref')
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
    cdnas = []
    if (not ngs.utils.all_exists(cdna_path, t2g_path)) or overwrite:
        for fasta_path, gtf_path in zip(fasta_paths, gtf_paths):
            logger.info(f'Preparing {fasta_path}, {gtf_path}')
            # Parse GTF for gene and transcripts
            gene_infos, transcript_infos = ngs.gtf.genes_and_transcripts_from_gtf(
                gtf_path, use_version=True
            )

            # Split
            cdna_temp_path = get_temporary_filename(temp_dir)
            logger.info(
                f'Splitting genome {fasta_path} into cDNA at {cdna_temp_path}'
            )
            cdna_temp_path = ngs.fasta.split_genomic_fasta_to_cdna(
                fasta_path, cdna_temp_path, gene_infos, transcript_infos
            )
            cdnas.append(cdna_temp_path)

        logger.info(f'Concatenating {len(cdnas)} cDNAs to {cdna_path}')
        cdna_path = concatenate_files(*cdnas, out_path=cdna_path)
        results.update({'cdna_fasta': cdna_path})
    else:
        logger.info(
            f'Skipping cDNA FASTA generation because {cdna_path} already exists. Use --overwrite flag to overwrite'
        )

    if not glob.glob(f'{index_path}*') or overwrite:
        t2g_result = create_t2g_from_fasta(cdna_path, t2g_path)
        results.update(t2g_result)

        if k and k != 31:
            logger.warning(
                f'Using provided k-mer length {k} instead of optimal length 31'
            )
        index_result = split_and_index(
            cdna_path, index_path, n=n, k=k or 31, temp_dir=temp_dir
        ) if n > 1 else kallisto_index(
            cdna_path, index_path, k=k or 31
        )
        results.update(index_result)
    else:
        logger.info(
            'Skipping kallisto index because {} already exists. Use the --overwrite flag to overwrite.'
            .format(index_path)
        )

    return results


@logger.namespaced('ref_kite')
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


@logger.namespaced('ref_lamanno')
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
    if (not ngs.utils.all_exists(cdna_path, intron_path, t2g_path,
                                 cdna_t2c_path, intron_t2c_path)) or overwrite:
        for fasta_path, gtf_path in zip(fasta_paths, gtf_paths):
            logger.info(f'Preparing {fasta_path}, {gtf_path}')
            # Parse GTF for gene and transcripts
            gene_infos, transcript_infos = ngs.gtf.genes_and_transcripts_from_gtf(
                gtf_path, use_version=True
            )

            # Split cDNA
            cdna_temp_path = get_temporary_filename(temp_dir)
            logger.info(
                f'Splitting genome {fasta_path} into cDNA at {cdna_temp_path}'
            )
            cdna_temp_path = ngs.fasta.split_genomic_fasta_to_cdna(
                fasta_path, cdna_temp_path, gene_infos, transcript_infos
            )
            cdnas.append(cdna_temp_path)

            # cDNA t2c
            cdna_t2c_temp_path = get_temporary_filename(temp_dir)
            logger.info(
                f'Creating cDNA transcripts-to-capture at {cdna_t2c_temp_path}'
            )
            cdna_t2c_result = create_t2c(cdna_temp_path, cdna_t2c_temp_path)
            cdna_t2cs.append(cdna_t2c_result['t2c'])

            # Split intron
            intron_temp_path = get_temporary_filename(temp_dir)
            logger.info(f'Splitting genome into introns at {intron_temp_path}')
            intron_temp_path = ngs.fasta.split_genomic_fasta_to_intron(
                fasta_path,
                intron_temp_path,
                gene_infos,
                transcript_infos,
                flank=flank if flank is not None else k -
                1 if k is not None else 30
            )
            introns.append(intron_temp_path)

            # intron t2c
            intron_t2c_temp_path = get_temporary_filename(temp_dir)
            logger.info(
                f'Creating intron transcripts-to-capture at {intron_t2c_temp_path}'
            )
            intron_t2c_result = create_t2c(
                intron_temp_path, intron_t2c_temp_path
            )
            intron_t2cs.append(intron_t2c_result['t2c'])

        # Concatenate
        logger.info(f'Concatenating {len(cdnas)} cDNA FASTAs to {cdna_path}')
        cdna_path = concatenate_files(*cdnas, out_path=cdna_path)
        logger.info(
            f'Concatenating {len(cdna_t2cs)} cDNA transcripts-to-captures to {cdna_t2c_path}'
        )
        cdna_t2c_path = concatenate_files(*cdna_t2cs, out_path=cdna_t2c_path)
        logger.info(
            f'Concatenating {len(introns)} intron FASTAs to {intron_path}'
        )
        intron_path = concatenate_files(*introns, out_path=intron_path)
        logger.info(
            f'Concatenating {len(intron_t2cs)} intron transcripts-to-captures to {intron_t2c_path}'
        )
        intron_t2c_path = concatenate_files(
            *intron_t2cs, out_path=intron_t2c_path
        )
        results.update({
            'cdna_fasta': cdna_path,
            'cdna_t2c': cdna_t2c_path,
            'intron_fasta': intron_path,
            'intron_t2c': intron_t2c_path
        })

    else:
        logger.info(
            'Skipping cDNA and intron FASTA generation because files already exist. Use --overwrite flag to overwrite'
        )

    if not glob.glob(f'{index_path}*') or overwrite:
        # Concatenate cDNA and intron fastas to generate T2G and build index
        combined_path = get_temporary_filename(temp_dir)
        logger.info(f'Concatenating cDNA and intron FASTAs to {combined_path}')
        combined_path = concatenate_files(
            cdna_path, intron_path, out_path=combined_path
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

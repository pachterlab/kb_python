import logging
import os
import tarfile
from urllib.request import urlretrieve

from .config import (
    get_kallisto_binary_path,
    INDEX_FILENAME,
    REFERENCES_MAPPING,
    T2G_FILENAME,
)
from .constants import (
    COMBINED_FILENAME,
    SORTED_FASTA_FILENAME,
    SORTED_GTF_FILENAME,
)
from .fasta import (
    FASTA,
    generate_cdna_fasta,
    generate_intron_fasta,
)
from .gtf import GTF
from .utils import (
    concatenate_files,
    open_as_text,
    run_executable,
)

logger = logging.getLogger(__name__)


def sort_gtf(gtf_path, out_path):
    """Sorts a GTF file based on its chromosome, start position, line number.

    :param gtf_path: path to GTF file
    :type gtf_path: str

    :return: path to sorted GTF file
    :rtype: str
    """
    logger.info('Sorting {}'.format(gtf_path))
    gtf = GTF(gtf_path)
    gtf.sort(out_path)
    return out_path


def sort_fasta(fasta_path, out_path):
    """Sorts a FASTA file based on its header.

    :param fasta_path: path to FASTA file
    :type fasta_path: str

    :return: path to sorted FASTA file
    :rtype: str
    """
    logger.info('Sorting {}'.format(fasta_path))
    fasta = FASTA(fasta_path)
    fasta.sort(out_path)
    return out_path


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
                gene_name = entry['group'].get('gene_name', '')
                f.write('{}\t{}\t{}\n'.format(transcript, gene, gene_name))

                if intron:
                    f.write(
                        '{}\t{}\t{}\n'.format(
                            transcript + '-I', gene, gene_name
                        )
                    )

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
        for sequence_id, _ in fasta.entries():
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
    logger.info('Indexing to {}'.format(index_path))
    command = [
        get_kallisto_binary_path(), 'index', '-i', index_path, '-k', k,
        fasta_path
    ]
    run_executable(command)
    return {'index': index_path}


def download_reference(
        choice, index_path, t2g_path, temp_dir='tmp', overwrite=False
):
    """Downloads a provided reference file from a static url.

    The configuration for provided references is in `config.py`.

    :param choice: reference key
    :type choice: str
    :param index_path: path to output kallisto index
    :type index_path: str
    :param t2g_path: path to output transcript-to-gene mapping
    :type t2g_path: str
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional
    :param overwrite: overwrite an existing index file, defaults to `False`
    :type overwrite: bool, optional

    :return: dictionary containing paths to generated file(s)
    :rtype: dict
    """
    results = {}
    if not os.path.exists(index_path) or overwrite:
        reference = REFERENCES_MAPPING[choice]
        url = reference.url
        path = os.path.join(temp_dir, os.path.basename(url))
        logging.info(
            'Downloading files for {} from {} to {}'.format(choice, url, path)
        )
        local_path, headers = urlretrieve(url, path)

        logging.info('Extracting files from {}'.format(local_path))
        with tarfile.open(local_path, 'r:gz') as f:
            f.extractall(temp_dir)

        os.rename(os.path.join(temp_dir, INDEX_FILENAME), index_path)
        results.update({'index': index_path})
        os.rename(os.path.join(temp_dir, T2G_FILENAME), t2g_path)
        results.update({'t2g': t2g_path})
    else:
        logger.info(
            'Skipping download because {} already exists. Use the --overwrite flag to overwrite.'
            .format(index_path)
        )
    return results


def ref(
        fasta_path,
        gtf_path,
        cdna_path,
        index_path,
        t2g_path,
        temp_dir='tmp',
        overwrite=False
):
    """Generates files necessary to generate count matrices for single-cell RNA-seq.

    :param fasta_path: path to genomic FASTA file
    :type fasta_path: str
    :param gtf_path: path to GTF file
    :type gtf_path: str
    :param cdna_path: path to generate the cDNA FASTA file
    :type cdna_path: str
    :param t2g_path: path to output transcript-to-gene mapping
    :type t2g_path: str
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional
    :param overwrite: overwrite an existing index file, defaults to `False`
    :type overwrite: bool, optional

    :return: dictionary containing paths to generated file(s)
    :rtype: dict
    """
    results = {}
    t2g_result = create_t2g_from_gtf(gtf_path, t2g_path)
    results.update(t2g_result)
    if not os.path.exists(index_path) or overwrite:
        sorted_fasta_path = sort_fasta(
            fasta_path, os.path.join(temp_dir, SORTED_FASTA_FILENAME)
        )
        sorted_gtf_path = sort_gtf(
            gtf_path, os.path.join(temp_dir, SORTED_GTF_FILENAME)
        )
        logger.info('Splitting genome into cDNA at {}'.format(cdna_path))
        cdna_fasta_path = generate_cdna_fasta(
            sorted_fasta_path, sorted_gtf_path, cdna_path
        )

        index_result = kallisto_index(cdna_fasta_path, index_path)
        results.update(index_result)
    else:
        logger.info(
            'Skipping kallisto index because {} already exists. Use the --overwrite flag to overwrite.'
            .format(index_path)
        )
    return results


def ref_lamanno(
        fasta_path,
        gtf_path,
        cdna_path,
        intron_path,
        index_path,
        t2g_path,
        cdna_t2c_path,
        intron_t2c_path,
        temp_dir='tmp',
        overwrite=False,
):
    """Generates files necessary to generate RNA velocity matrices for single-cell RNA-seq.

    :param fasta_path: path to genomic FASTA file
    :type fasta_path: str
    :param gtf_path: path to GTF file
    :type gtf_path: str
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
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional
    :param overwrite: overwrite an existing index file, defaults to `False`
    :type overwrite: bool, optional

    :return: dictionary containing paths to generated file(s)
    :rtype: dict
    """
    results = {}
    t2g_result = create_t2g_from_gtf(gtf_path, t2g_path, intron=True)
    results.update(t2g_result)
    if not os.path.exists(index_path) or overwrite:
        sorted_fasta_path = sort_fasta(
            fasta_path, os.path.join(temp_dir, SORTED_FASTA_FILENAME)
        )
        sorted_gtf_path = sort_gtf(
            gtf_path, os.path.join(temp_dir, SORTED_GTF_FILENAME)
        )
        logger.info('Splitting genome into cDNA at {}'.format(cdna_path))
        cdna_fasta_path = generate_cdna_fasta(
            sorted_fasta_path, sorted_gtf_path, cdna_path
        )
        results.update({'cdna_fasta': cdna_fasta_path})
        logger.info(
            'Creating cDNA transcripts-to-capture at {}'.format(cdna_t2c_path)
        )
        cdna_t2c_result = create_t2c(cdna_fasta_path, cdna_t2c_path)
        results.update({'cdna_t2c': cdna_t2c_result['t2c']})
        logger.info('Splitting genome into introns at {}'.format(intron_path))
        intron_fasta_path = generate_intron_fasta(
            sorted_fasta_path, sorted_gtf_path, intron_path
        )
        results.update({'intron_fasta': intron_fasta_path})
        logger.info(
            'Creating intron transcripts-to-capture at {}'.
            format(cdna_t2c_path)
        )
        intron_t2c_result = create_t2c(intron_fasta_path, intron_t2c_path)
        results.update({'intron_t2c': intron_t2c_result['t2c']})
        logger.info('Concatenating cDNA and intron FASTAs')
        combined_path = concatenate_files(
            cdna_fasta_path,
            intron_fasta_path,
            out_path=os.path.join(temp_dir, COMBINED_FILENAME),
            temp_dir=temp_dir
        )
        index_result = kallisto_index(combined_path, index_path)
        results.update(index_result)
    else:
        logger.info(
            'Skipping kallisto index because {} already exists. Use the --overwrite flag to overwrite.'
            .format(index_path)
        )

    return results

import concurrent.futures
import logging
import os
import re
import shutil
import subprocess as sp
import threading
import time
from urllib.request import urlretrieve

import anndata
import ngs_tools as ngs
import pandas as pd
import scipy.io
from scipy import sparse
from tqdm import tqdm

from .config import (
    get_bustools_binary_path,
    get_kallisto_binary_path,
    PLATFORM,
    TECHNOLOGIES_MAPPING,
    UnsupportedOSException,
)
from .dry import dryable
from .dry import utils as dry_utils
from .logging import logger
from .stats import STATS

TECHNOLOGY_PARSER = re.compile(r'^(?P<name>\S+)')
VERSION_PARSER = re.compile(r'^\S*? ([0-9]+).([0-9]+).([0-9]+)')

# These functions have been moved as of 0.26.1 to the ngs_tools library but are
# imported from this file in other places. For now, let's keep these here.
# TODO: remove these
open_as_text = ngs.utils.open_as_text
decompress_gzip = ngs.utils.decompress_gzip
compress_gzip = ngs.utils.compress_gzip
concatenate_files = ngs.utils.concatenate_files_as_text
download_file = ngs.utils.download_file
get_temporary_filename = dryable(dry_utils.get_temporary_filename)(
    ngs.utils.mkstemp
)


class NotImplementedException(Exception):
    pass


class UnmetDependencyException(Exception):
    pass


class TqdmLoggingHandler(logging.Handler):
    """Custom logging handler so that logging does not affect progress bars.
    """

    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception:
            self.handleError(record)


def update_filename(filename, code):
    """Update the provided path with the specified code.

    For instance, if the `path` is 'output.bus' and `code` is `s` (for sort),
    this function returns `output.s.bus`.

    :param filename: filename (NOT path)
    :type filename: str
    :param code: code to append to filename
    :type code: str

    :return: path updated with provided code
    :rtype: str
    """
    name, extension = os.path.splitext(filename)
    return f'{name}.{code}{extension}'


@dryable(dry_utils.make_directory)
def make_directory(path):
    """Quietly make the specified directory (and any subdirectories).

    This function is a wrapper around os.makedirs. It is used so that
    the appropriate mkdir command can be printed for dry runs.

    :param path: path to directory to make
    :type path: str
    """
    os.makedirs(path, exist_ok=True)


@dryable(dry_utils.remove_directory)
def remove_directory(path):
    """Quietly make the specified directory (and any subdirectories).

    This function is a wrapper around shutil.rmtree. It is used so that
    the appropriate rm command can be printed for dry runs.

    :param path: path to directory to remove
    :type path: str
    """
    shutil.rmtree(path, ignore_errors=True)


@dryable(dry_utils.run_executable)
def run_executable(
    command,
    stdin=None,
    stdout=sp.PIPE,
    stderr=sp.PIPE,
    wait=True,
    stream=True,
    quiet=False,
    returncode=0,
    alias=True,
    record=True,
):
    """Execute a single shell command.

    :param command: a list representing a single shell command
    :type command: list
    :param stdin: object to pass into the `stdin` argument for `subprocess.Popen`,
                  defaults to `None`
    :type stdin: stream, optional
    :param stdout: object to pass into the `stdout` argument for `subprocess.Popen`,
                  defaults to `subprocess.PIPE`
    :type stdout: stream, optional
    :param stderr: object to pass into the `stderr` argument for `subprocess.Popen`,
                  defaults to `subprocess.PIPE`
    :type stderr: stream, optional
    :param wait: whether to wait until the command has finished, defaults to `True`
    :type wait: bool, optional
    :param stream: whether to stream the output to the command line, defaults to `True`
    :type stream: bool, optional
    :param quiet: whether to not display anything to the command line and not check the return code,
                  defaults to `False`
    :type quiet: bool, optional
    :param returncode: the return code expected if the command runs as intended,
                       defaults to `0`
    :type returncode: int, optional
    :param alias: whether to use the basename of the first element of `command`,
                  defaults to `True`
    :type alias: bool, optional
    :param record: whether to record the call statistics, defaults to `True`
    :type record: bool, optional

    :return: the spawned process
    :rtype: subprocess.Process
    """
    command = [str(c) for c in command]
    c = command.copy()
    if alias:
        c[0] = os.path.basename(c[0])
    if not quiet:
        logger.debug(' '.join(c))
    if not wait and record:
        STATS.command(c)
    start = time.time()
    p = sp.Popen(
        command,
        stdin=stdin,
        stdout=stdout,
        stderr=stderr,
        universal_newlines=wait,
        bufsize=1 if wait else -1,
    )

    # Wait if desired.
    if wait:
        out = []
        while p.poll() is None:
            if stream and not quiet:
                for line in p.stdout:
                    out.append(line.strip())
                    logger.debug(line.strip())
                for line in p.stderr:
                    out.append(line.strip())
                    logger.debug(line.strip())
        if record:
            STATS.command(c, runtime=time.time() - start)

        if not quiet and p.returncode != returncode:
            logger.error('\n'.join(out))
            raise sp.CalledProcessError(p.returncode, ' '.join(command))

    return p


def get_kallisto_version():
    """Get the provided Kallisto version.

    This function parses the help text by executing the included Kallisto binary.

    :return: tuple of major, minor, patch versions
    :rtype: tuple
    """
    p = run_executable([get_kallisto_binary_path()],
                       quiet=True,
                       returncode=1,
                       record=False)
    match = VERSION_PARSER.match(p.stdout.read())
    return tuple(int(ver) for ver in match.groups()) if match else None


def get_bustools_version():
    """Get the provided Bustools version.

    This function parses the help text by executing the included Bustools binary.

    :return: tuple of major, minor, patch versions
    :rtype: tuple
    """
    p = run_executable([get_bustools_binary_path()],
                       quiet=True,
                       returncode=1,
                       record=False)
    match = VERSION_PARSER.match(p.stdout.read())
    return tuple(int(ver) for ver in match.groups()) if match else None


def parse_technologies(lines):
    """Parse a list of strings into a list of supported technologies.

    This function parses the technologies printed by running `kallisto bus --list`.

    :param lines: the output of `kallisto bus --list` split into lines
    :type lines: list

    :return: list of technologies
    :rtype: list
    """
    parsing = False
    technologies = set()
    for line in lines:
        if line.startswith('-'):
            parsing = True
            continue

        if parsing:
            if line.isspace():
                break
            match = TECHNOLOGY_PARSER.match(line)
            if match:
                technologies.add(match['name'])
    return technologies


def get_supported_technologies():
    """Runs 'kallisto bus --list' to fetch a list of supported technologies.

    :return: list of technologies
    :rtype: list
    """
    p = run_executable([get_kallisto_binary_path(), 'bus', '--list'],
                       quiet=True,
                       returncode=1,
                       record=False)
    return parse_technologies(p.stdout)


def whitelist_provided(technology):
    """Determine whether or not the whitelist for a technology is provided.

    :param technology: the name of the technology
    :type technology: str

    :return: whether the whitelist is provided
    :rtype: bool
    """
    upper = technology.upper()
    return upper in TECHNOLOGIES_MAPPING and TECHNOLOGIES_MAPPING[
        upper].chemistry.has_whitelist


@dryable(dry_utils.move_file)
def move_file(source, destination):
    """Move a file from source to destination, overwriting the file if the
    destination exists.

    :param source: path to source file
    :type source: str
    :param destination: path to destination
    :type destination: str

    :return: path to moved file
    :rtype: str
    """
    shutil.move(source, destination)
    return destination


@dryable(dry_utils.copy_whitelist)
def copy_whitelist(technology, out_dir):
    """Copies provided whitelist for specified technology.

    :param technology: the name of the technology
    :type technology: str
    :param out_dir: directory to put the whitelist
    :type out_dir: str

    :return: path to whitelist
    :rtype: str
    """
    technology = TECHNOLOGIES_MAPPING[technology.upper()]
    archive_path = technology.chemistry.whitelist_path
    whitelist_path = os.path.join(
        out_dir,
        os.path.splitext(os.path.basename(archive_path))[0]
    )
    with open_as_text(archive_path, 'r') as f, open(whitelist_path, 'w') as out:
        out.write(f.read())
    return whitelist_path


@dryable(dry_utils.copy_map)
def copy_map(technology, out_dir):
    """Copies provided feature-to-cell barcode mapping for the speified technology.

    :param technology: the name of the technology
    :type technology: str
    :param out_dir: directory to put the map
    :type out_dir: str

    :return: path to map
    :rtype: str
    """
    technology = TECHNOLOGIES_MAPPING[technology.upper()]
    archive_path = technology.chemistry.whitelist_path
    map_path = os.path.join(
        out_dir,
        os.path.splitext(os.path.basename(archive_path))[0]
    )
    with open_as_text(archive_path, 'r') as f, open(map_path, 'w') as out:
        out.write(f.read())
    return map_path


@dryable(dry_utils.stream_file)
def stream_file(url, path):
    """Creates a FIFO file to use for piping remote files into processes.

    This function spawns a new thread to download the remote file into a FIFO
    file object. FIFO file objects are only supported on unix systems.

    :param url: url to the file
    :type url: str
    :param path: path to place FIFO file
    :type path: str

    :raises UnsupportedOSException: if the OS is Windows

    :return: path to FIFO file
    :rtype: str
    """
    # Windows does not support FIFO files.
    if PLATFORM == 'windows':
        raise UnsupportedOSException((
            'Windows does not support piping remote files.'
            'Please download the file manually.'
        ))
    else:
        logger.info('Piping {} to {}'.format(url, path))
        os.mkfifo(path)
        t = threading.Thread(target=urlretrieve, args=(url, path), daemon=True)
        t.start()
    return path


def read_t2g(t2g_path):
    """Given a transcript-to-gene mapping path, read it into a dictionary.
    The first column is always assumed to tbe the transcript IDs.

    :param t2g_path: path to t2g
    :type t2g_path: str

    :return: dictionary containing transcript IDs as keys and all other columns
             as a tuple as values
    :rtype: dict
    """
    t2g = {}
    with open_as_text(t2g_path, 'r') as f:
        for line in f:
            if line.isspace():
                continue
            split = line.strip().split('\t')
            transcript = split[0]
            other = tuple(split[1:])
            if transcript in t2g:
                logger.warning(
                    f'Found duplicate entries for {transcript} in {t2g_path}. '
                    'Earlier entries will be ignored.'
                )
            t2g[transcript] = other
    return t2g


def import_tcc_matrix_as_anndata(
    matrix_path, barcodes_path, ec_path, txnames_path, threads=8
):
    """Import a TCC matrix as an Anndata object.

    :param matrix_path: path to the matrix ec file
    :type matrix_path: str
    :param barcodes_path: path to the barcodes txt file
    :type barcodes_path: str
    :param genes_path: path to the ec txt file
    :type genes_path: str
    :param txnames_path: path to transcripts.txt generated by `kallisto bus`
    :type txnames_path: str

    :return: a new Anndata object
    :rtype: anndata.Anndata
    """
    df_barcodes = pd.read_csv(
        barcodes_path, index_col=0, header=None, names=['barcode']
    )
    df_ec = pd.read_csv(
        ec_path,
        index_col=0,
        header=None,
        names=['ec', 'transcripts'],
        sep='\t',
        dtype=str
    )
    df_ec.index = df_ec.index.astype(str)  # To prevent logging from anndata
    with open(txnames_path, 'r') as f:
        transcripts = [
            line.strip() for line in f.readlines() if not line.strip().isspace()
        ]

    ts = list(df_ec.transcripts)
    get_transcript_ids = lambda ts, transcripts: [
        [transcripts[int(i)] for i in t.split(',')] for t in ts
    ]
    futures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        chunk = int(len(ts) / threads) + 1
        for i in range(threads):
            future = executor.submit(
                get_transcript_ids, ts[i * chunk:(i + 1) * chunk], transcripts
            )
            futures.append(future)
    transcript_ids = []
    for future in futures:
        transcript_ids += future.result()
    df_ec['transcript_ids'] = transcript_ids
    df_ec.drop('transcripts', axis=1, inplace=True)
    return anndata.AnnData(
        X=scipy.io.mmread(matrix_path).tocsr(), obs=df_barcodes, var=df_ec
    )


def import_matrix_as_anndata(
    matrix_path, barcodes_path, genes_path, t2g_path=None, name='gene'
):
    """Import a matrix as an Anndata object.

    :param matrix_path: path to the matrix ec file
    :type matrix_path: str
    :param barcodes_path: path to the barcodes txt file
    :type barcodes_path: str
    :param genes_path: path to the genes txt file
    :type genes_path: str
    :param t2g_path: path to transcript-to-gene mapping. If this is provided,
                     the third column of the mapping is appended to the
                     anndata var, defaults to `None`
    :type t2g_path: str, optional
    :param name: name of the columns, defaults to "gene"
    :type name: str, optional

    :return: a new Anndata object
    :rtype: anndata.Anndata
    """
    df_barcodes = pd.read_csv(
        barcodes_path, index_col=0, header=None, names=['barcode']
    )
    df_genes = pd.read_csv(
        genes_path, header=None, index_col=0, names=[f'{name}_id'], sep='\t'
    )
    df_genes.index = df_genes.index.astype(
        str
    )  # To prevent logging from anndata
    mtx = scipy.io.mmread(matrix_path)

    # If any of the genes are duplicated, collapse them by summing
    if any(df_genes.index.duplicated()):
        logger.debug(
            f'Deduplicating genes found in {genes_path} by adding duplicates'
        )
        mtx = mtx.tocsc()
        gene_indices = {}
        for i, gene in enumerate(df_genes.index):
            gene_indices.setdefault(gene, []).append(i)
        genes = []
        deduplicated_mtx = sparse.lil_matrix(
            (len(df_barcodes), len(gene_indices))
        )
        for i, gene in enumerate(sorted(gene_indices.keys())):
            genes.append(gene)
            indices = gene_indices[gene]
            deduplicated_mtx[:, i] = mtx[:, indices].sum(axis=1)
        df_genes = pd.DataFrame(index=pd.Series(genes, name=f'{name}_id'))
        mtx = deduplicated_mtx

    t2g = read_t2g(t2g_path) if t2g_path else {}
    id_to_name = {}
    for transcript, attributes in t2g.items():
        if len(attributes) > 1:
            id_to_name[attributes[0]] = attributes[1]
    gene_names = [id_to_name.get(i, '') for i in df_genes.index]
    if any(bool(g) for g in gene_names):
        df_genes[f'{name}_name'] = gene_names

    return anndata.AnnData(X=mtx.tocsr(), obs=df_barcodes, var=df_genes)


def overlay_anndatas(adata_spliced, adata_unspliced):
    """'Overlays' anndata objects by taking the intersection of the obs and var
    of each anndata.

    :param adata_spliced: an Anndata object
    :type adata_spliced: anndata.Anndata
    :param adata_unspliced: an Anndata object
    :type adata_unspliced: anndata.Anndata

    :return: a new Anndata object
    :rtype: anndata.Anndata
    """
    obs_idx = adata_spliced.obs.index.intersection(adata_unspliced.obs.index)
    var_idx = adata_spliced.var.index.intersection(adata_unspliced.var.index)
    spliced_intersection = adata_spliced[obs_idx][:, var_idx]
    unspliced_intersection = adata_unspliced[obs_idx][:, var_idx]

    df_obs = unspliced_intersection.obs
    df_var = unspliced_intersection.var
    return anndata.AnnData(
        X=spliced_intersection.X,
        layers={
            'spliced': spliced_intersection.X,
            'unspliced': unspliced_intersection.X
        },
        obs=df_obs,
        var=df_var
    )


def sum_anndatas(adata_spliced, adata_unspliced):
    """Sum the counts in two anndata objects by taking the intersection of
    both matrices and adding the values together.

    :param adata_spliced: an Anndata object
    :type adata_spliced: anndata.Anndata
    :param adata_unspliced: an Anndata object
    :type adata_unspliced: anndata.Anndata

    :return: a new Anndata object
    :rtype: anndata.Anndata
    """
    obs_idx = adata_spliced.obs.index.intersection(adata_unspliced.obs.index)
    var_idx = adata_spliced.var.index.intersection(adata_unspliced.var.index)
    spliced_intersection = adata_spliced[obs_idx][:, var_idx]
    unspliced_intersection = adata_unspliced[obs_idx][:, var_idx]

    df_obs = unspliced_intersection.obs
    df_var = unspliced_intersection.var
    return anndata.AnnData(
        X=spliced_intersection.X + unspliced_intersection.X,
        obs=df_obs,
        var=df_var
    )

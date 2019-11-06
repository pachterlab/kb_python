import gzip
import logging
import os
import re
import shutil
import subprocess as sp
import threading
import time
from urllib.request import urlretrieve

import anndata
import pandas as pd
import scipy.io

from .config import (
    get_bustools_binary_path,
    get_kallisto_binary_path,
    PACKAGE_PATH,
    PLATFORM,
    TECHNOLOGIES_MAPPING,
    WHITELIST_DIR,
    UnsupportedOSException,
)

logger = logging.getLogger(__name__)

TECHNOLOGY_PARSER = re.compile(r'^(?P<name>\S+)')
VERSION_PARSER = re.compile(r'^\S*? ([0-9]+).([0-9]+).([0-9]+)')


class NotImplementedException(Exception):
    pass


class UnmetDependencyException(Exception):
    pass


def open_as_text(path, mode):
    """Open a textfile or gzip file in text mode.

    :param path: path to textfile or gzip
    :type path: str
    :param mode: mode to open the file, either `w` for write or `r` for read
    :type mode: str

    :return: file object
    :rtype: file object
    """
    return gzip.open(path, mode +
                     't') if path.endswith('.gz') else open(path, mode)


def decompress_gzip(gzip_path, out_path):
    """Decompress a gzip file to provided file path.

    :param gzip_path: path to gzip file
    :type gzip_path: str
    :param out_path: path to decompressed file
    :type out_path: str

    :return: path to decompressed file
    :rtype: str
    """
    with gzip.open(gzip_path, 'rb') as f, open(out_path, 'wb') as out:
        shutil.copyfileobj(f, out)
    return out_path


def compress_gzip(file_path, out_path):
    """Compress a file into gzip.

    :param file_path: path to file
    :type file_path: str
    :param out_dir: path to compressed file
    :type out_dir: str

    :return: path to compressed file
    :rtype: str
    """
    with open(file_path, 'rb') as f, gzip.open(out_path, 'wb') as out:
        shutil.copyfileobj(f, out)
    return out_path


def run_executable(
        command,
        stdin=None,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        wait=True,
        stream=True,
        quiet=False,
        returncode=0
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

    :return: the spawned process
    :rtype: subprocess.Process
    """
    command = [str(c) for c in command]
    if not quiet:
        logger.debug(' '.join(command))
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
        while p.poll() is None:
            if stream and not quiet:
                for line in p.stdout:
                    logger.debug(line.strip())
                for line in p.stderr:
                    logger.debug(line.strip())
            else:
                time.sleep(1)

        if not quiet and p.returncode != returncode:
            raise sp.CalledProcessError(p.returncode, ' '.join(command))

    return p


def run_chain(*commands, stdin=None, stdout=sp.PIPE, wait=True, stream=False):
    """Execute multiple shell commands by piping the output into inputs.

    :param commands: lists of shell commands
    :type commands: list
    :param stdin: object to pass into the `stdin` argument for `subprocess.Popen`,
                  defaults to `None`
    :type stdin: stream, optional
    :param stdout: object to pass into the `stdout` argument for `subprocess.Popen`,
                  defaults to `subprocess.PIPE`
    :type stdout: stream, optional
    :param wait: whether to wait until the command has finished, defaults to `True`
    :type wait: bool, optional
    :param stream: whether to stream the output to the command line, defaults to `True`
    :type stream: bool, optional

    :return: list of spawned subprocesses
    :rtype: list
    """
    assert len(commands) > 1
    processes = []
    for command in commands:
        _stdin = stdin
        _stdout = stdout
        _wait = wait
        _stream = stream
        if processes:
            _stdin = processes[-1].stdout
        if len(processes) != len(commands) - 1:
            _stdout = sp.PIPE
            _wait = False
            _stream = False
        p = run_executable(
            command, stdin=_stdin, stdout=_stdout, wait=_wait, stream=_stream
        )
        processes.append(p)

        if _stdin:
            _stdin.close()

    if wait:
        for p in processes:
            while p.poll() is None:
                time.sleep(1)
            if p.returncode != 0:
                raise sp.CalledProcessError(p.returncode, ' '.join(command))
    return processes


def get_kallisto_version():
    """Get the provided Kallisto version.

    This function parses the help text by executing the included Kallisto binary.

    :return: tuple of major, minor, patch versions
    :rtype: tuple
    """
    p = run_executable([get_kallisto_binary_path()], quiet=True, returncode=1)
    match = VERSION_PARSER.match(p.stdout.read())
    return tuple(int(ver) for ver in match.groups()) if match else None


def get_bustools_version():
    """Get the provided Bustools version.

    This function parses the help text by executing the included Bustools binary.

    :return: tuple of major, minor, patch versions
    :rtype: tuple
    """
    p = run_executable([get_bustools_binary_path()], quiet=True, returncode=1)
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
                       returncode=1)
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
        upper].whitelist_archive


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
    archive_path = os.path.join(
        PACKAGE_PATH, WHITELIST_DIR, technology.whitelist_archive
    )
    whitelist_path = os.path.join(
        out_dir,
        os.path.splitext(technology.whitelist_archive)[0]
    )
    with open_as_text(archive_path, 'r') as f, open(whitelist_path, 'w') as out:
        out.write(f.read())
    return whitelist_path


def concatenate_files(*paths, out_path, temp_dir='tmp'):
    """Concatenates an arbitrary number of files into one TEXT file.

    Only supports text and gzip files.

    :param paths: an arbitrary number of paths to files
    :type paths: str
    :param out_path: path to place concatenated file
    :type out_path: str
    :param temp_dir: temporary directory, defaults to `tmp`
    :type temp_dir: str, optional

    :return: path to concatenated file
    :rtype: str
    """
    with open(out_path, 'w') as out:
        for path in paths:
            with open_as_text(path, 'r') as f:
                for line in f:
                    if not line.isspace():
                        out.write(line.strip() + '\n')

    return out_path


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


def import_matrix_as_anndata(matrix_path, barcodes_path, genes_path):
    """Import a matrix as an Anndata object.

    :param matrix_path: path to the matrix ec file
    :type matrix_path: str
    :param barcodes_path: path to the barcodes txt file
    :type barcodes_path: str
    :param genes_path: path to the genes txt file
    :type genes_path: str

    :return: a new Anndata object
    :rtype: anndata.Anndata
    """
    df_barcodes = pd.read_csv(
        barcodes_path, index_col=0, header=None, names=['barcode']
    )
    df_genes = pd.read_csv(
        genes_path, header=None, index_col=0, names=['gene_id'], sep='\t'
    )
    return anndata.AnnData(
        X=scipy.io.mmread(matrix_path).tocsr(), obs=df_barcodes, var=df_genes
    )


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
    idx = adata_spliced.obs.index.intersection(adata_unspliced.obs.index)
    spliced_intersection = adata_spliced[idx]
    unspliced_intersection = adata_unspliced[idx]
    spliced_unspliced = spliced_intersection.copy()
    spliced_unspliced.layers['spliced'] = spliced_intersection.X
    spliced_unspliced.layers['unspliced'] = unspliced_intersection.X
    return spliced_unspliced


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
    idx = adata_spliced.obs.index.intersection(adata_unspliced.obs.index)
    spliced_intersection = adata_spliced[idx]
    unspliced_intersection = adata_unspliced[idx]
    spliced_unspliced = spliced_intersection.copy()
    spliced_unspliced.X = spliced_intersection.X + unspliced_intersection.X
    return spliced_unspliced

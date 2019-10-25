import gzip
import logging
import os
import re
import subprocess as sp
import tarfile
import time

import anndata
import pandas as pd
import scipy.io

from .config import (
    get_bustools_binary_path,
    get_kallisto_binary_path,
    PACKAGE_PATH,
    TECHNOLOGIES_MAPPING,
    WHITELIST_DIR,
)

logger = logging.getLogger(__name__)

TECHNOLOGY_PARSER = re.compile(r'^(?P<name>\S+)')
VERSION_PARSER = re.compile(r'^\S*? ([0-9]+).([0-9]+).([0-9]+)')


class NotImplementedException(Exception):
    pass


class UnmetDependencyException(Exception):
    pass


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
    """Run a single shell command and wait for it to terminate.
    """
    command = [str(c) for c in command]
    if not quiet:
        logger.info(' '.join(command))
    p = sp.Popen(
        command,
        stdin=stdin,
        stdout=stdout,
        stderr=stderr,
        universal_newlines=wait
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
    """Chain multiple commands by piping outputs to inputs.
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
    p = run_executable([get_kallisto_binary_path()], quiet=True, returncode=1)
    match = VERSION_PARSER.match(p.stdout.read())
    return tuple(int(ver) for ver in match.groups()) if match else None


def get_bustools_version():
    p = run_executable([get_bustools_binary_path()], quiet=True, returncode=1)
    match = VERSION_PARSER.match(p.stdout.read())
    return tuple(int(ver) for ver in match.groups()) if match else None


def parse_technologies(lines):
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
    """
    p = run_executable([get_kallisto_binary_path(), 'bus', '--list'],
                       quiet=True,
                       returncode=1)
    return parse_technologies(p.stdout)


def whitelist_provided(technology):
    upper = technology.upper()
    return upper in TECHNOLOGIES_MAPPING and TECHNOLOGIES_MAPPING[
        upper].whitelist_archive


def copy_whitelist(technology, out_dir):
    """Copies provided whitelist barcodes for specified technology.
    """
    if not TECHNOLOGIES_MAPPING[technology.upper()]:
        raise NotImplementedException(
            'whitelist for {} is not yet provided by kb_python'.
            format(technology)
        )

    technology = TECHNOLOGIES_MAPPING[technology.upper()]
    archive_path = os.path.join(
        PACKAGE_PATH, WHITELIST_DIR, technology.whitelist_archive
    )
    with tarfile.open(archive_path, 'r:gz') as f:
        f.extract(technology.whitelist_filename, path=out_dir)
    return os.path.join(out_dir, technology.whitelist_filename)


def concatenate_files(*paths, out_path, temp_dir='tmp'):
    """Concatenates an arbitrary number of files into one TEXT file.

    Only supports text and gzip files.
    """
    with open(out_path, 'w') as out:
        for path in paths:
            with gzip.open(path, 'rt') if path.endswith('.gz') else open(
                    path, 'r') as f:
                for line in f:
                    if not line.isspace():
                        out.write(line.strip() + '\n')

    return out_path


def import_matrix_as_anndata(matrix_path, barcodes_path, genes_path):
    df_barcodes = pd.read_csv(
        barcodes_path, index_col=0, header=None, names=['barcode']
    )
    df_genes = pd.read_csv(
        genes_path, header=None, index_col=0, names=['gene_id'], sep='\t'
    )
    df_genes.index = df_genes.index.str.split('.').str[
        0]  # slice off version number
    return anndata.AnnData(
        X=scipy.io.mmread(matrix_path).tocsr(), obs=df_barcodes, var=df_genes
    )


def overlay_anndatas(*adatas):
    pass


def convert_matrix_to_loom(matrix_path, barcodes_path, genes_path, out_path):
    adata = import_matrix_as_anndata(matrix_path, barcodes_path, genes_path)
    adata.write_loom(out_path)
    return out_path


def convert_matrix_to_h5ad(matrix_path, barcodes_path, genes_path, out_path):
    adata = import_matrix_as_anndata(matrix_path, barcodes_path, genes_path)
    adata.write(out_path)
    return out_path

import logging
import os
import re
import subprocess as sp
import tarfile
import time

import anndata
import pandas as pd
import scipy.io

from .config import TECHNOLOGIES_MAPPING, WHITELIST_DIR
from .constants import MINIMUM_REQUIREMENTS

logger = logging.getLogger(__name__)

TECHNOLOGY_PARSER = re.compile(r'^(?P<name>\S+)')
VERSION_PARSERS = {
    requirement:
    re.compile(r'^{} ([0-9]+).([0-9]+).([0-9]+)'.format(requirement))
    for requirement in MINIMUM_REQUIREMENTS
}


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
        stream=False,
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
            if stream:
                for line in p.stdout:
                    logger.info(line.strip())
                for line in p.stdout:
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


def get_version(requirement):
    p = run_executable([requirement], quiet=True, returncode=1)
    match = VERSION_PARSERS[requirement].match(p.stdout.read())
    return tuple(int(ver) for ver in match.groups()) if match else None


def check_dependencies():
    """Checks if executable dependencies have been met.
    """
    for requirement, minimum_version in MINIMUM_REQUIREMENTS.items():
        version = get_version(requirement)
        if version < minimum_version:
            raise UnmetDependencyException(
                '{} version {} is less than the minimum requirement {}'.format(
                    requirement, version, minimum_version
                )
            )


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
    p = run_executable(['kallisto', 'bus', '--list'], quiet=True, returncode=1)
    return parse_technologies(p.stdout)


def copy_whitelist(technology):
    """Copies provided whitelist barcodes for specified technology.
    """
    if not TECHNOLOGIES_MAPPING[technology.upper()]:
        raise NotImplementedException(
            'whitelist for {} is not yet provided by kb_python'.
            format(technology)
        )

    technology = TECHNOLOGIES_MAPPING[technology.upper()]
    archive_path = os.path.join(
        os.path.dirname(__file__), WHITELIST_DIR, technology.whitelist_archive
    )
    whitelist_filename = technology.whitelist_filename
    with tarfile.open(archive_path, 'r:gz') as f:
        f.extract(whitelist_filename)
    return whitelist_filename


def create_transcript_list(gtf_path, use_name=True, use_version=True):
    r = {}
    with open(gtf_path, 'r') as f:
        for line in f:
            if len(line) == 0 or line[0] == '#':
                continue
            l = line.strip().split('\t')  # noqa
            if l[2] == 'transcript':
                info = l[8]
                d = {}
                for x in info.split('; '):
                    x = x.strip()
                    p = x.find(' ')
                    if p == -1:
                        continue
                    k = x[:p]
                    p = x.find('"', p)
                    p2 = x.find('"', p + 1)
                    v = x[p + 1:p2]
                    d[k] = v

                if 'transcript_id' not in d or 'gene_id' not in d:
                    continue

                tid = d['transcript_id'].split(".")[0]
                gid = d['gene_id'].split(".")[0]
                if use_version:
                    if 'transcript_version' not in d or 'gene_version' not in d:
                        continue

                    tid += '.' + d['transcript_version']
                    gid += '.' + d['gene_version']
                gname = None
                if use_name:
                    if 'gene_name' not in d:
                        continue
                    gname = d['gene_name']

                if tid in r:
                    continue

                r[tid] = (gid, gname)
    return r


def import_matrix_as_anndata(matrix_path, barcodes_path, genes_path):
    df_barcodes = pd.read_csv(
        barcodes_path, index_col=0, header=None, names=['barcode']
    )
    df_genes = pd.read_csv(
        genes_path, header=None, index_col=0, names=['ensembl_id'], sep='\t'
    )
    df_genes.index = df_genes.index.str.slice(0, 18)  # slice off version number
    return anndata.AnnData(
        X=scipy.io.mmread(matrix_path).tocsr(), obs=df_barcodes, var=df_genes
    )


def convert_matrix_to_loom(matrix_path, barcodes_path, genes_path, out_path):
    adata = import_matrix_as_anndata(matrix_path, barcodes_path, genes_path)
    adata.write_loom(out_path)
    return out_path


def convert_matrix_to_h5ad(matrix_path, barcodes_path, genes_path, out_path):
    adata = import_matrix_as_anndata(matrix_path, barcodes_path, genes_path)
    adata.write(out_path)
    return out_path

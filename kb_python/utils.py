import concurrent.futures
import functools
import os
import queue
import re
import shutil
import subprocess as sp
import threading
import time
from typing import Callable, Dict, List, Optional, Set, Tuple, Union
from urllib.request import urlretrieve

import anndata
import ngs_tools as ngs
import pandas as pd
import scipy.io
from scipy import sparse

from .config import (
    get_bustools_binary_path,
    get_kallisto_binary_path,
    PLATFORM,
    TECHNOLOGIES_MAPPING,
    UnsupportedOSError,
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


def update_filename(filename: str, code: str) -> str:
    """Update the provided path with the specified code.

    For instance, if the `path` is 'output.bus' and `code` is `s` (for sort),
    this function returns `output.s.bus`.

    Args:
        filename: filename (NOT path)
        code: code to append to filename

    Returns:
        Path updated with provided code
    """
    name, extension = os.path.splitext(filename)
    return f'{name}.{code}{extension}'


@dryable(dry_utils.make_directory)
def make_directory(path: str):
    """Quietly make the specified directory (and any subdirectories).

    This function is a wrapper around os.makedirs. It is used so that
    the appropriate mkdir command can be printed for dry runs.

    Args:
        path: Path to directory to make
    """
    os.makedirs(path, exist_ok=True)


@dryable(dry_utils.remove_directory)
def remove_directory(path: str):
    """Quietly make the specified directory (and any subdirectories).

    This function is a wrapper around shutil.rmtree. It is used so that
    the appropriate rm command can be printed for dry runs.

    Args:
        path: Path to directory to remove
    """
    shutil.rmtree(path, ignore_errors=True)


@dryable(dry_utils.run_executable)
def run_executable(
    command: List[str],
    stdin: Optional[int] = None,
    stdout: int = sp.PIPE,
    stderr: int = sp.PIPE,
    wait: bool = True,
    stream: bool = True,
    quiet: bool = False,
    returncode: int = 0,
    alias: bool = True,
    record: bool = True,
) -> Union[Tuple[sp.Popen, str, str], sp.Popen]:
    """Execute a single shell command.

    Args:
        command: A list representing a single shell command
        stdin: Object to pass into the `stdin` argument for `subprocess.Popen`,
            defaults to `None`
        stdout: Object to pass into the `stdout` argument for `subprocess.Popen`,
            defaults to `subprocess.PIPE`
        stderr: Object to pass into the `stderr` argument for `subprocess.Popen`,
            defaults to `subprocess.PIPE`
        wait: Whether to wait until the command has finished, defaults to `True`
        stream: Whether to stream the output to the command line, defaults to `True`
        quiet: Whether to not display anything to the command line and not check the return code,
            defaults to `False`
        returncode: The return code expected if the command runs as intended,
            defaults to `0`
        alias: Whether to use the basename of the first element of `command`,
            defaults to `True`
        record: Whether to record the call statistics, defaults to `True`

    Returns:
        (the spawned process, list of strings printed to stdout,
            list of strings printed to stderr) if `wait=True`.
            Otherwise, the spawned process
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

    # Helper function to read from a pipe and put the output to a queue.
    def reader(pipe, qu, stop_event, name):
        while not stop_event.is_set():
            for _line in pipe:
                line = _line.strip()
                qu.put((name, line))

    # Wait if desired.
    if wait:
        stdout = ''
        stderr = ''
        out = []
        out_queue = queue.Queue()
        stop_event = threading.Event()
        stdout_reader = threading.Thread(
            target=reader,
            args=(p.stdout, out_queue, stop_event, 'stdout'),
            daemon=True
        )
        stderr_reader = threading.Thread(
            target=reader,
            args=(p.stderr, out_queue, stop_event, 'stderr'),
            daemon=True
        )
        stdout_reader.start()
        stderr_reader.start()

        while p.poll() is None:
            while not out_queue.empty():
                name, line = out_queue.get()
                if stream and not quiet:
                    logger.debug(line)
                out.append(line)
                if name == 'stdout':
                    stdout += f'{line}\n'
                elif name == 'stderr':
                    stderr += f'{line}\n'
            else:
                time.sleep(0.1)

        # Stop readers & flush queue
        stop_event.set()
        time.sleep(1)
        while not out_queue.empty():
            name, line = out_queue.get()
            if stream and not quiet:
                logger.debug(line)
            out.append(line)
            if name == 'stdout':
                stdout += f'{line}\n'
            elif name == 'stderr':
                stderr += f'{line}\n'
        if record:
            STATS.command(c, runtime=time.time() - start)

        if not quiet and p.returncode != returncode:
            logger.error('\n'.join(out))
            raise sp.CalledProcessError(p.returncode, ' '.join(command))

    return (p, stdout, stderr) if wait else p


def get_kallisto_version() -> Optional[Tuple[int, int, int]]:
    """Get the provided Kallisto version.

    This function parses the help text by executing the included Kallisto binary.

    Returns:
        Major, minor, patch versions
    """
    p, stdout, stderr = run_executable([get_kallisto_binary_path()],
                                       quiet=True,
                                       returncode=1,
                                       record=False)
    match = VERSION_PARSER.match(stdout)
    return tuple(int(ver) for ver in match.groups()) if match else None


def get_bustools_version() -> Optional[Tuple[int, int, int]]:
    """Get the provided Bustools version.

    This function parses the help text by executing the included Bustools binary.

    Returns:
        Major, minor, patch versions
    """
    p, stdout, stderr = run_executable([get_bustools_binary_path()],
                                       quiet=True,
                                       returncode=1,
                                       record=False)
    match = VERSION_PARSER.match(stdout)
    return tuple(int(ver) for ver in match.groups()) if match else None


def parse_technologies(lines: List[str]) -> Set[str]:
    """Parse a list of strings into a list of supported technologies.

    This function parses the technologies printed by running `kallisto bus --list`.

    Args:
        lines: The output of `kallisto bus --list` split into lines

    Returns:
        Set of technologies
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


def get_supported_technologies() -> Set[str]:
    """Runs 'kallisto bus --list' to fetch a list of supported technologies.

    Returns:
        Set of technologies
    """
    p, stdout, stderr = run_executable([
        get_kallisto_binary_path(), 'bus', '--list'
    ],
                                       quiet=True,
                                       returncode=1,
                                       record=False)
    return parse_technologies(stdout)


def whitelist_provided(technology: str) -> bool:
    """Determine whether or not the whitelist for a technology is provided.

    Args:
        technology: The name of the technology

    Returns:
        Whether the whitelist is provided
    """
    upper = technology.upper()
    return upper in TECHNOLOGIES_MAPPING and TECHNOLOGIES_MAPPING[
        upper].chemistry.has_whitelist


@dryable(dry_utils.move_file)
def move_file(source: str, destination: str) -> str:
    """Move a file from source to destination, overwriting the file if the
    destination exists.

    Args:
        source: Path to source file
        destination: Path to destination

    Returns:
        Path to moved file
    """
    shutil.move(source, destination)
    return destination


@dryable(dry_utils.copy_whitelist)
def copy_whitelist(technology: str, out_dir: str) -> str:
    """Copies provided whitelist for specified technology.

    Args:
        technology: The name of the technology
        out_dir: Directory to put the whitelist

    Returns:
        Path to whitelist
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
def copy_map(technology: str, out_dir: str) -> str:
    """Copies provided feature-to-cell barcode mapping for the speified technology.

    Args:
        technology: The name of the technology
        out_dir: Directory to put the map

    Returns:
        Path to map
    """
    technology = TECHNOLOGIES_MAPPING[technology.upper()]
    archive_path = technology.chemistry.feature_map_path
    map_path = os.path.join(
        out_dir,
        os.path.splitext(os.path.basename(archive_path))[0]
    )
    with open_as_text(archive_path, 'r') as f, open(map_path, 'w') as out:
        out.write(f.read())
    return map_path


@dryable(dry_utils.stream_file)
def stream_file(url: str, path: str) -> str:
    """Creates a FIFO file to use for piping remote files into processes.

    This function spawns a new thread to download the remote file into a FIFO
    file object. FIFO file objects are only supported on unix systems.

    Args:
        url: Url to the file
        path: Path to place FIFO file

    Returns:
        Path to FIFO file

    Raises:
        UnsupportedOSError: If the OS is Windows
    """
    # Windows does not support FIFO files.
    if PLATFORM == 'windows':
        raise UnsupportedOSError((
            'Windows does not support piping remote files.'
            'Please download the file manually.'
        ))
    else:
        logger.info('Piping {} to {}'.format(url, path))
        os.mkfifo(path)
        t = threading.Thread(target=urlretrieve, args=(url, path), daemon=True)
        t.start()
    return path


def read_t2g(t2g_path: str) -> Dict[str, Tuple[str, ...]]:
    """Given a transcript-to-gene mapping path, read it into a dictionary.
    The first column is always assumed to tbe the transcript IDs.

    Args:
        t2g_path: Path to t2g

    Returns:
        Dictionary containing transcript IDs as keys and all other columns
            as a tuple as values
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


def collapse_anndata(
    adata: anndata.AnnData, by: Optional[str] = None
) -> anndata.AnnData:
    """Collapse the given Anndata by summing duplicate rows. The `by` argument
    specifies which column to use. If not provided, the index is used.

    Note:
        This function also collapses any existing layers. Additionally, the
        returned AnnData will have the values used to collapse as the index.

    Args:
        adata: The Anndata to collapse
        by: The column to collapse by. If not provided, the index is used. When
            this column contains missing values (i.e. nan or None), these
            columns are removed.

    Returns:
        A new collapsed Anndata object. All matrices are sparse, regardless of
        whether or not they were in the input Anndata.
    """
    var = adata.var
    if by is not None:
        var = var.set_index(by)
    na_mask = var.index.isna()
    adata = adata[:, ~na_mask].copy()
    adata.var = var[~na_mask]

    if not any(adata.var.index.duplicated()):
        return adata

    var_indices = {}
    for i, index in enumerate(adata.var.index):
        var_indices.setdefault(index, []).append(i)

    # Convert all original matrices to csc for fast column operations
    X = sparse.csc_matrix(adata.X)
    layers = {
        layer: sparse.csc_matrix(adata.layers[layer])
        for layer in adata.layers
    }
    new_index = []
    # lil_matrix is efficient for row-by-row construction
    new_X = sparse.lil_matrix((len(var_indices), adata.shape[0]))
    new_layers = {layer: new_X.copy() for layer in adata.layers}
    for i, (index, indices) in enumerate(var_indices.items()):
        new_index.append(index)
        new_X[i] = X[:, indices].sum(axis=1).flatten()
        for layer in layers.keys():
            new_layers[layer][i] = layers[layer][:,
                                                 indices].sum(axis=1).flatten()

    return anndata.AnnData(
        X=new_X.T.tocsr(),
        layers={layer: new_layers[layer].T.tocsr()
                for layer in new_layers},
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=pd.Series(new_index, name=adata.var.index.name)),
    )


def import_tcc_matrix_as_anndata(
    matrix_path: str,
    barcodes_path: str,
    ec_path: str,
    txnames_path: str,
    threads: int = 8
) -> anndata.AnnData:
    """Import a TCC matrix as an Anndata object.

    Args:
        matrix_path: Path to the matrix ec file
        barcodes_path: Path to the barcodes txt file
        genes_path: Path to the ec txt file
        txnames_path: Path to transcripts.txt generated by `kallisto bus`

    Returns:
        A new Anndata object
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
        ';'.join(transcripts[int(i)] for i in t.split(',')) for t in ts
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
    df_ec['transcript_ids'] = pd.Categorical(transcript_ids)
    df_ec.drop('transcripts', axis=1, inplace=True)
    return anndata.AnnData(
        X=scipy.io.mmread(matrix_path).tocsr(), obs=df_barcodes, var=df_ec
    )


def import_matrix_as_anndata(
    matrix_path: str,
    barcodes_path: str,
    genes_path: str,
    t2g_path: Optional[str] = None,
    name: str = 'gene',
    by_name: bool = False,
) -> anndata.AnnData:
    """Import a matrix as an Anndata object.

    Args:
        matrix_path: Path to the matrix ec file
        barcodes_path: Path to the barcodes txt file
        genes_path: Path to the genes txt file
        t2g_path: Path to transcript-to-gene mapping. If this is provided,
            the third column of the mapping is appended to the anndata var,
            defaults to `None`
        name: Name of the columns, defaults to "gene"
        by_name: Aggregate counts by name instead of ID. `t2g_path` must be
            provided and contain names.

    Returns:
        A new Anndata object
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
    adata = collapse_anndata(
        anndata.AnnData(X=mtx.tocsr(), obs=df_barcodes, var=df_genes)
    )

    name_column = f'{name}_name'
    if t2g_path:
        t2g = read_t2g(t2g_path)
        id_to_name = {}
        for transcript, attributes in t2g.items():
            if len(attributes) > 1:
                id_to_name[attributes[0]] = attributes[1]
        gene_names = [id_to_name.get(i, '') for i in adata.var.index]
        if any(bool(g) for g in gene_names):
            adata.var[name_column] = pd.Categorical(gene_names)

    return (
        collapse_anndata(adata, by=name_column)
        if name_column in adata.var.columns and by_name else adata
    )


def overlay_anndatas(
    adata_spliced: anndata.AnnData, adata_unspliced: anndata.AnnData
) -> anndata.AnnData:
    """'Overlays' anndata objects by taking the intersection of the obs and var
    of each anndata.

    Note:
        Matrices generated by kallisto | bustools always contain all genes,
        even if they have zero counts. Therefore, taking the intersection
        is not entirely necessary but is done as a sanity check.

    Args:
        adata_spliced: An Anndata object
        adata_unspliced: An Anndata object

    Returns:
        A new Anndata object
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


def sum_anndatas(
    adata_spliced: anndata.AnnData, adata_unspliced: anndata.AnnData
) -> anndata.AnnData:
    """Sum the counts in two anndata objects by taking the intersection of
    both matrices and adding the values together.

    Note:
        Matrices generated by kallisto | bustools always contain all genes,
        even if they have zero counts. Therefore, taking the intersection
        is not entirely necessary but is done as a sanity check.

    Args:
        adata_spliced: An Anndata object
        adata_unspliced: An Anndata object

    Returns:
        A new Anndata object
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


def restore_cwd(func: Callable) -> Callable:
    """Function decorator to decorate functions that change the current working
    directory. When such a function is decorated with this function, the
    current working directory is restored to its previous state when the
    function exits.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        old_cwd = os.path.abspath(os.getcwd())
        try:
            return func(*args, **kwargs)
        finally:
            os.chdir(old_cwd)

    return wrapper

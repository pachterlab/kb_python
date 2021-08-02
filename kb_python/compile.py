import os
import shutil
import tempfile
from urllib.parse import urljoin, urlparse

import requests

from .config import (
    BUSTOOLS_RELEASES_URL,
    BUSTOOLS_TARBALL_URL,
    COMPILED_DIR,
    KALLISTO_RELEASES_URL,
    KALLISTO_TARBALL_URL,
    PLATFORM,
)
from .logging import logger
from .utils import download_file, restore_cwd, run_executable


def get_latest_github_release_tag(releases_url):
    """Get the tag name of the latest GitHub release, given a url to the
    releases API.

    :param releases_url: url to the releases API
    :type releases_url: str

    :return: tag name
    :rtype: str
    """
    response = requests.get(releases_url)
    response.raise_for_status()
    return response.json()[0]['tag_name']


def get_filename_from_url(url):
    response = requests.get(url)
    response.raise_for_status()
    disposition = response.headers.get('content-disposition')
    if disposition:
        for split in disposition.split(';'):
            split = split.strip()
            if split.startswith('filename'):
                return split[split.index('=') + 1:].strip('\"\'')
    else:
        return os.path.basename(urlparse(url).path)


def get_kallisto_url():
    """Get the tarball url of the latest kallisto release.

    :return: tarball url
    :rtype: str
    """
    tag = get_latest_github_release_tag(KALLISTO_RELEASES_URL)
    return urljoin(KALLISTO_TARBALL_URL, tag)


def get_bustools_url():
    """Get the tarball url of the latest bustools release.

    :return: tarball url
    :rtype: str
    """
    tag = get_latest_github_release_tag(BUSTOOLS_RELEASES_URL)
    return urljoin(BUSTOOLS_TARBALL_URL, tag)


def find_git_root(path):
    """Find the root directory of a git repo by walking.

    :param path: path to start the search
    :type path: str

    :return: path to root of git repo
    :rtype: str
    """
    for root, dirs, files in os.walk(path):
        if '.gitignore' in files:
            return root
    raise Exception('Unable to find git root.')


@restore_cwd
def compile_kallisto(source_dir, binary_path, cmake_arguments=None):
    """Compile `kallisto` from source.

    :param source_dir: path to directory containing root of kallisto git repo
    :type source_dir: str
    :param binary_path: path to place compiled binary
    :type binary_path: str
    :param cmake_arguments: additional arguments to pass to the cmake command
    :type cmake_arguments: str, optional

    :return: path to compiled binary
    :rtype: str
    """
    source_dir = os.path.abspath(source_dir)
    binary_path = os.path.abspath(binary_path)
    os.makedirs(os.path.dirname(binary_path), exist_ok=True)

    logger.info(
        f'Compiling `kallisto` binary from source at {source_dir} to {binary_path}. '
        'This requires `autoheader`, `autoconf`, `cmake` and `make` to be executable '
        'from the command-line, as well as zlib development headers. '
        'See https://pachterlab.github.io/kallisto/source for more information.'
    )
    os.chdir(source_dir)
    shutil.copyfile(
        'license.txt',
        os.path.join(os.path.dirname(binary_path), 'license.txt')
    )

    os.chdir(os.path.join('ext', 'htslib'))
    run_executable(['autoheader'])
    run_executable(['autoconf'])
    os.chdir(os.path.join('..', '..'))
    os.makedirs('build', exist_ok=True)
    os.chdir('build')
    cmake_command = ['cmake', '..']
    if cmake_arguments:
        cmake_command.append(cmake_arguments)
    run_executable(cmake_command)
    run_executable(['make'])
    os.makedirs(os.path.dirname(binary_path), exist_ok=True)
    shutil.copy2(
        os.path.join(
            'src', 'kallisto.exe' if PLATFORM == 'windows' else 'kallisto'
        ), binary_path
    )
    return binary_path


@restore_cwd
def compile_bustools(source_dir, binary_path, cmake_arguments=None):
    """Compile `bustools` from source.

    :param source_dir: path to directory containing root of bustools git repo
    :type source_dir: str
    :param binary_path: path to place compiled binary
    :type binary_path: str
    :param cmake_arguments: additional arguments to pass to the cmake command
    :type cmake_arguments: str, optional

    :return: path to compiled binary
    :rtype: str
    """
    source_dir = os.path.abspath(source_dir)
    binary_path = os.path.abspath(binary_path)
    os.makedirs(os.path.dirname(binary_path), exist_ok=True)

    logger.info(
        f'Compiling `bustools` binary from source {source_dir} to {binary_path}. '
        'This requires `cmake` and `make` to be executable from the command-line. '
        'See https://bustools.github.io/source for more information.'
    )
    os.chdir(source_dir)
    shutil.copyfile(
        'LICENSE', os.path.join(os.path.dirname(binary_path), 'LICENSE')
    )

    os.makedirs('build', exist_ok=True)
    os.chdir('build')
    cmake_command = ['cmake', '..']
    if cmake_arguments:
        cmake_command.append(cmake_arguments)
    run_executable(cmake_command)
    run_executable(['make'])
    shutil.copy2(
        os.path.join(
            'src', 'bustools.exe' if PLATFORM == 'windows' else 'bustools'
        ), binary_path
    )
    return binary_path


@logger.namespaced('compile')
def compile(
    target,
    out_dir=None,
    cmake_arguments=None,
    url=None,
    overwrite=False,
    temp_dir='tmp',
):
    """Compile `kallisto` and/or `bustools` binaries by downloading and compiling
    a source archive.

    :param target: which binary to compile. May be one of `kallisto`, `bustools`
        or `all`
    :type target: str
    :param out_dir: path to output directory, defaults to `None`
    :type out_dir: str, optional
    :param cmake_arguments: additional arguments to pass to the cmake command
    :type cmake_arguments: str, optional
    :param url: download the source archive from this url instead, defaults to
        `None`
    :type url: str, optional
    :param overwrite: overwrite any existing results, defaults to `False`
    :type overwrite: bool, optional
    :param temp_dir: path to temporary directory, defaults to `tmp`
    :type temp_dir: str, optional

    :return: dictionary of results
    :rtype: dict
    """
    results = {}
    if target in ('kallisto', 'all'):
        binary_path = os.path.join(
            out_dir or os.path.join(COMPILED_DIR, 'kallisto'), 'kallisto'
        )
        if os.path.exists(binary_path) and not overwrite:
            raise Exception(
                f'Compiled binary already exists at {binary_path}. '
                'Use `--overwrite` to overwrite.'
            )

        _url = url or get_kallisto_url()
        logger.info(f'Downloading kallisto source from {_url}')
        archive_path = download_file(
            _url, os.path.join(temp_dir, get_filename_from_url(_url))
        )
        source_dir = tempfile.mkdtemp(dir=temp_dir)
        shutil.unpack_archive(archive_path, source_dir)
        source_dir = find_git_root(source_dir)
        binary_path = compile_kallisto(
            source_dir, binary_path, cmake_arguments=cmake_arguments
        )
        results['kallisto'] = binary_path
    if target in ('bustools', 'all'):
        binary_path = os.path.join(
            out_dir or os.path.join(COMPILED_DIR, 'bustools'), 'bustools'
        )
        if os.path.exists(binary_path) and not overwrite:
            raise Exception(
                f'Compiled binary already exists at {binary_path}. '
                'Use `--overwrite` to overwrite.'
            )

        _url = url or get_bustools_url()
        logger.info(f'Downloading bustools source from {_url}')
        archive_path = download_file(
            _url, os.path.join(temp_dir, get_filename_from_url(_url))
        )
        source_dir = tempfile.mkdtemp(dir=temp_dir)
        shutil.unpack_archive(archive_path, source_dir)
        source_dir = find_git_root(source_dir)
        binary_path = compile_bustools(
            source_dir, binary_path, cmake_arguments=cmake_arguments
        )
        results['bustools'] = binary_path
    return results

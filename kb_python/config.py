import os
import platform
from collections import namedtuple

PACKAGE_PATH = os.path.dirname(__file__)
PLATFORM = platform.system().lower()
BINS_DIR = 'bins'

TEMP_DIR = 'tmp'

# Technology to file position mapping
Technology = namedtuple(
    'Technology', [
        'name', 'description', 'nfiles', 'reads_file', 'umi_positions',
        'barcode_positions', 'whitelist_archive'
    ]
)
WHITELIST_DIR = 'whitelists'
TECHNOLOGIES = [
    Technology(
        '10XV1',
        '10x version 1',
        3,
        0,
        [(1, 0, 0)],
        [(2, 0, 0)],
        '10xv1_whitelist.txt.gz',
    ),
    Technology(
        '10XV2', '10x version 2', 2, 1, [(0, 16, 26)], [(0, 0, 16)],
        '10xv2_whitelist.txt.gz'
    ),
    Technology(
        '10XV3', '10x version 3', 2, 1, [(0, 16, 28)], [(0, 0, 16)],
        '10xv3_whitelist.txt.gz'
    ),
    Technology('CELSEQ', 'CEL-Seq', 2, 1, [(0, 8, 12)], [(0, 0, 8)], None),
    Technology(
        'CELSEQ2', 'CEL-SEQ version 2', 2, 1, [(0, 0, 6)], [(0, 6, 12)], None
    ),
    Technology('DROPSEQ', 'DropSeq', 2, 1, [(0, 12, 20)], [(0, 0, 12)], None),
    Technology(
        'INDROPS',
        'inDrops version 3',
        2,
        1,
        [(0, 42, 48)],
        [(0, 0, 11), (0, 30, 38)],
        'inDropsv3_whitelist.txt.gz',
    ),
    Technology('SCRUBSEQ', 'SCRB-Seq', 2, 1, [(0, 6, 16)], [(0, 0, 6)], None),
    Technology(
        'SURECELL', 'SureCell for ddSEQ', 2, 1, [(0, 51, 59)], [(0, 0, 6),
                                                                (0, 21, 27),
                                                                (0, 42, 48)],
        None
    ),
]
TECHNOLOGIES_MAPPING = {t.name: t for t in TECHNOLOGIES}

# Supported pre-built indices
Reference = namedtuple(
    'Reference', ['name', 'url', 'lamanno_url', 'velocity_url']
)
REFERENCES = [
    Reference(
        'human',
        'https://caltech.box.com/shared/static/v1nm7lpnqz5syh8dyzdk2zs8bglncfib.gz',
        None, None
    ),
    Reference(
        'mouse',
        'https://caltech.box.com/shared/static/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz',
        None, None
    )
]
REFERENCES_MAPPING = {r.name: r for r in REFERENCES}
# File names that are in the tar.gz file.
INDEX_FILENAME = 'transcriptome.idx'
T2G_FILENAME = 'transcripts_to_genes.txt'


class UnsupportedOSException(Exception):
    pass


def get_kallisto_binary_path():
    """Get the path to the platform-dependent Kallisto binary included with
    the installation.

    :return: path to the binary
    :rtype: str
    """
    bin_filename = 'kallisto.exe' if PLATFORM == 'windows' else 'kallisto'
    path = os.path.join(
        PACKAGE_PATH, BINS_DIR, PLATFORM, 'kallisto', bin_filename
    )
    if not os.path.exists(path):
        raise UnsupportedOSException(
            'This operating system ({}) is not supported.'.format(PLATFORM)
        )
    return path


def get_bustools_binary_path():
    """Get the path to the platform-dependent Bustools binary included with
    the installation.

    :return: path to the binary
    :rtype: str
    """
    bin_filename = 'bustools.exe' if PLATFORM == 'windows' else 'bustools'
    path = os.path.join(
        PACKAGE_PATH, BINS_DIR, PLATFORM, 'bustools', bin_filename
    )
    if not os.path.exists(path):
        raise UnsupportedOSException(
            'This operating system ({}) is not supported.'.format(PLATFORM)
        )
    return path

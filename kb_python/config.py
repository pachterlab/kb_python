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
        'name', 'description', 'nfiles', 'seq_pos', 'umi_pos',
        'whitelist_archive', 'whitelist_filename'
    ]
)
WHITELIST_DIR = 'whitelists'
TECHNOLOGIES = [
    Technology(
        '10XV1', '10x version 1', 3, 0, 1, '10xv1_whitelist.tar.gz',
        '10xv1_whitelist.txt'
    ),
    Technology(
        '10XV2', '10x version 2', 2, 1, 0, '10xv2_whitelist.tar.gz',
        '10xv2_whitelist.txt'
    ),
    Technology(
        '10XV3', '10x version 3', 2, 1, 0, '10xv3_whitelist.tar.gz',
        '10xv3_whitelist.txt'
    ),
    Technology('CELSEQ', 'CEL-Seq', 2, 1, 0, None, None),
    Technology('CELSEQ2', 'CEL-SEQ version 2', 2, 1, 0, None, None),
    Technology('DROPSEQ', 'DropSeq', 2, 1, 0, None, None),
    Technology('INDROPS', 'inDrops', 2, 1, 0, None, None),
    Technology('SCRUBSEQ', 'SCRB-Seq', 2, 1, 0, None, None),
    Technology('SURECELL', 'SureCell for ddSEQ', 2, 1, 0, None, None),
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
    bin_filename = 'bustools.exe' if PLATFORM == 'windows' else 'bustools'
    path = os.path.join(
        PACKAGE_PATH, BINS_DIR, PLATFORM, 'bustools', bin_filename
    )
    if not os.path.exists(path):
        raise UnsupportedOSException(
            'This operating system ({}) is not supported.'.format(PLATFORM)
        )
    return path

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
        'name', 'nfiles', 'seq_pos', 'umi_pos', 'whitelist_archive',
        'whitelist_filename'
    ]
)
WHITELIST_DIR = 'whitelists'
TECHNOLOGIES = [
    Technology(
        '10XV1', 3, 0, 1, '10xv1_whitelist.tar.gz', '10xv1_whitelist.txt'
    ),
    Technology(
        '10XV2', 2, 1, 0, '10xv2_whitelist.tar.gz', '10xv2_whitelist.txt'
    ),
    Technology(
        '10XV3', 2, 1, 0, '10xv3_whitelist.tar.gz', '10xv3_whitelist.txt'
    ),
    Technology('CELSEQ', 2, 1, 0, None, None),
    Technology('CELSEQ2', 2, 1, 0, None, None),
    Technology('DROPSEQ', 2, 1, 0, None, None),
    Technology('INDROPS', 2, 1, 0, None, None),
    Technology('SCRUBSEQ', 2, 1, 0, None, None),
    Technology('SURECELL', 2, 1, 0, None, None),
]
TECHNOLOGIES_MAPPING = {t.name: t for t in TECHNOLOGIES}


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

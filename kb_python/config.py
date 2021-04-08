import os
import platform
import shutil
from collections import namedtuple

PACKAGE_PATH = os.path.abspath(os.path.dirname(__file__))
PLATFORM = platform.system().lower()
BINS_DIR = 'bins'

TEMP_DIR = 'tmp'
DRY = False
VALIDATE = True
CHUNK_SIZE = 1024 * 1024 * 4  # Download files in chunks of 4 Mb


def get_provided_kallisto_path():
    """Finds platform-dependent kallisto binary included with the installation.

    :return: path to the binary, `None` if not found
    :rtype: str
    """
    bin_filename = 'kallisto.exe' if PLATFORM == 'windows' else 'kallisto'
    path = os.path.join(
        PACKAGE_PATH, BINS_DIR, PLATFORM, 'kallisto', bin_filename
    )
    if not os.path.isfile(path):
        return None
    return path


def get_provided_bustools_path():
    """Finds platform-dependent bustools binary included with the installation.

    :return: path to the binary, `None` if not found
    :rtype: str
    """
    bin_filename = 'bustools.exe' if PLATFORM == 'windows' else 'bustools'
    path = os.path.join(
        PACKAGE_PATH, BINS_DIR, PLATFORM, 'bustools', bin_filename
    )
    if not os.path.isfile(path):
        return None
    return path


# Binary paths. These should hold the full path to the binaries that should
# be called throughout the execution of the program. Therefore, this
# usually needs to be set only once. Defaults to provided binaries.
KALLISTO_PATH = get_provided_kallisto_path()
BUSTOOLS_PATH = get_provided_bustools_path()

# Technology to file position mapping
Technology = namedtuple(
    'Technology', [
        'name', 'description', 'nfiles', 'reads_file', 'umi_positions',
        'barcode_positions', 'whitelist_archive', 'map_archive'
    ]
)
WHITELIST_DIR = 'whitelists'
MAP_DIR = 'maps'
TECHNOLOGIES = [
    Technology(
        '10XV1',
        '10x version 1',
        3,
        2,
        [(1, 0, 10)],
        [(0, 0, 14)],
        '10xv1_whitelist.txt.gz',
        None,
    ),
    Technology(
        '10XV2',
        '10x version 2',
        2,
        1,
        [(0, 16, 26)],
        [(0, 0, 16)],
        '10xv2_whitelist.txt.gz',
        None,
    ),
    Technology(
        '10XV3',
        '10x version 3',
        2,
        1,
        [(0, 16, 28)],
        [(0, 0, 16)],
        '10xv3_whitelist.txt.gz',
        '10xv3_feature_barcode_map.txt.gz',
    ),
    Technology(
        'CELSEQ',
        'CEL-Seq',
        2,
        1,
        [(0, 8, 12)],
        [(0, 0, 8)],
        None,
        None,
    ),
    Technology(
        'CELSEQ2',
        'CEL-SEQ version 2',
        2,
        1,
        [(0, 0, 6)],
        [(0, 6, 12)],
        None,
        None,
    ),
    Technology(
        'DROPSEQ',
        'DropSeq',
        2,
        1,
        [(0, 12, 20)],
        [(0, 0, 12)],
        None,
        None,
    ),
    Technology(
        'INDROPSV1',
        'inDrops version 1',
        2,
        1,
        [(0, 42, 48)],
        [(0, 0, 11), (0, 30, 38)],
        None,
        None,
    ),
    Technology(
        'INDROPSV2',
        'inDrops version 2',
        2,
        0,
        [(1, 42, 48)],
        [(1, 0, 11), (1, 30, 38)],
        None,
        None,
    ),
    Technology(
        'INDROPSV3',
        'inDrops version 3',
        3,
        2,
        [(1, 8, 14)],
        [(0, 0, 8), (1, 0, 8)],
        'inDropsv3_whitelist.txt.gz',
        None,
    ),
    Technology(
        'SCRUBSEQ',
        'SCRB-Seq',
        2,
        1,
        [(0, 6, 16)],
        [(0, 0, 6)],
        None,
        None,
    ),
    Technology(
        'SURECELL',
        'SureCell for ddSEQ',
        2,
        1,
        [(0, 51, 59)],
        [(0, 0, 6), (0, 21, 27), (0, 42, 48)],
        None,
        None,
    ),
    Technology(
        'SMARTSEQ', 'Smart-seq2', 2, '0, 1 (paired)', [], [], None, None
    )
]
TECHNOLOGIES_MAPPING = {t.name: t for t in TECHNOLOGIES}

# Supported pre-built indices
Reference = namedtuple('Reference', ['name', 'url', 'files'])
REFERENCES = [
    Reference(
        'human',
        'https://caltech.box.com/shared/static/v1nm7lpnqz5syh8dyzdk2zs8bglncfib.gz',
        {
            'i': 'transcriptome.idx',
            'g': 'transcripts_to_genes.txt'
        }
    ),
    Reference(
        'mouse',
        'https://caltech.box.com/shared/static/vcaz6cujop0xuapdmz0pplp3aoqc41si.gz',
        {
            'i': 'transcriptome.idx',
            'g': 'transcripts_to_genes.txt'
        }
    ),
    Reference(
        'linnarsson',
        'https://caltech.box.com/shared/static/kyf7ai5s8y2l0vycl5yxunrappvrf0yx.gz',
        {
            'i': 'gencode.v31.fragments.idx',
            'g': 'fragments2genes.txt',
            'c1': 'spliced_fragments.txt',
            'c2': 'unspliced_fragments.txt',
        }
    )
]
REFERENCES_MAPPING = {r.name: r for r in REFERENCES}


class UnsupportedOSException(Exception):
    pass


class NotExecutableException(Exception):
    pass


def get_kallisto_binary_path():
    """Dummy function that simply returns the current value of KALLISTO_PATH.
    """
    return KALLISTO_PATH


def get_bustools_binary_path():
    """Dummy function that simply returns the current value of BUSTOOLS_PATH.
    """
    return BUSTOOLS_PATH


def set_kallisto_binary_path(path):
    """Helper function to set the KALLISTO_PATH variable. Automatically finds the
    full path to the executable and sets that as KALLISTO_PATH.
    """
    global KALLISTO_PATH

    shutil_path = shutil.which(path)
    actual_path = None

    # First, check if it is an executable in the user's PATH
    if shutil_path:
        actual_path = os.path.abspath(shutil_path)
    elif os.path.isfile(path):
        actual_path = os.path.abspath(path)
    else:
        raise Exception(f'Unable to resolve path {path}')

    # Check that it is executable
    if not os.access(actual_path, os.X_OK):
        raise NotExecutableException(f'{actual_path} is not executable')

    KALLISTO_PATH = actual_path


def set_bustools_binary_path(path):
    """Helper function to set the BUSTOOLS_PATH variable. Automatically finds the
    full path to the executable and sets that as BUSTOOLS_PATH.
    """
    global BUSTOOLS_PATH

    shutil_path = shutil.which(path)
    actual_path = None

    # First, check if it is an executable in the user's PATH
    if shutil_path:
        actual_path = os.path.abspath(shutil_path)
    elif os.path.isfile(path):
        actual_path = os.path.abspath(path)
    else:
        raise Exception(f'Unable to resolve path {path}')

    # Check that it is executable
    if not os.access(actual_path, os.X_OK):
        raise NotExecutableException(f'{actual_path} is not executable')

    BUSTOOLS_PATH = actual_path


def set_dry():
    """Set this run to be a dry run.
    """
    global DRY
    DRY = True


def is_dry():
    """Return whether the current run is a dry run.

    :return: whether the current run is a dry run
    :rtype: bool
    """
    return DRY


def no_validate():
    """Turn off validation.
    """
    global VALIDATE
    VALIDATE = False


def is_validate():
    """Return whether validation is turned on.

    :return: whether validation is on
    :rtype: bool
    """
    return VALIDATE

import os
import tempfile

from ..config import (
    MAP_DIR,
    PACKAGE_PATH,
    PLATFORM,
    TECHNOLOGIES_MAPPING,
    WHITELIST_DIR,
    UnsupportedOSException,
)


def run_executable(command, quiet=False, alias=True, *args, **kwargs):
    """Dry version of `utils.run_executable`.
    """
    command = [str(c) for c in command]
    if not quiet:
        c = command.copy()
        if alias:
            c[0] = os.path.basename(c[0])
        print(' '.join(c))


def make_directory(path):
    """Dry version of `utils.make_directory`.
    """
    if PLATFORM == 'windows':
        print('md {}'.format(path))
    else:
        print('mkdir -p {}'.format(path))


def remove_directory(path):
    """Dry version of `utils.remove_directory`.
    """
    if PLATFORM == 'windows':
        print('rd /s /q "{}"'.format(path))
    else:
        print('rm -rf {}'.format(path))


def stream_file(url, path):
    """Dry version of `utils.stream_file`.
    """
    if PLATFORM == 'windows':
        raise UnsupportedOSException((
            'Windows does not support piping remote files.'
            'Please download the file manually.'
        ))
    else:
        print('mkfifo {}'.format(path))
        print('wget -bq {} -O {}'.format(url, path))


def move_file(source, destination):
    """Dry version of `utils.move_file`.
    """
    if PLATFORM == 'windows':
        print(f'move {source} {destination}')
    else:
        print(f'mv {source} {destination}')


def copy_whitelist(technology, out_dir):
    """Dry version of `utils.copy_whitelist`.
    """
    technology = TECHNOLOGIES_MAPPING[technology.upper()]
    archive_path = os.path.join(
        PACKAGE_PATH, WHITELIST_DIR, technology.whitelist_archive
    )
    whitelist_path = os.path.join(
        out_dir,
        os.path.splitext(technology.whitelist_archive)[0]
    )
    print('gzip -dc {} > {}'.format(archive_path, whitelist_path))
    return whitelist_path


def copy_map(technology, out_dir):
    """Dry version of `utils.copy_map`.
    """
    technology = TECHNOLOGIES_MAPPING[technology.upper()]
    archive_path = os.path.join(PACKAGE_PATH, MAP_DIR, technology.map_archive)
    map_path = os.path.join(
        out_dir,
        os.path.splitext(technology.map_archive)[0]
    )
    print('gzip -dc {} > {}'.format(archive_path, map_path))
    return map_path


def get_temporary_filename(temp_dir):
    """Dry version of `utils.get_temporary_filename`.
    """
    return os.path.join(
        temp_dir, f'{tempfile.gettempprefix()}{tempfile._get_candidate_names()}'
    )

import os

from ..config import PLATFORM, UnsupportedOSException


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

import functools
import os
import re

import scipy.io

from .config import get_bustools_binary_path, is_dry, is_validate
from .logging import logger
from .utils import run_executable

BUSTOOLS_INSPECT_PARSER = re.compile(r'^.*?(?P<count>[0-9]+)')


class FileVerificationFailed(Exception):
    pass


def validate_bus(path):
    """Verify if the provided BUS file is valid.

    A BUS file is considered valid when `bustools inspect` can read
    the file + it has > 0 BUS records.

    :param path: path to BUS file
    :type path: str

    :raises FileVerificationFailed: if the file failed verification
    :raises subprocess.CalledProcessError: if the bustools command failed
    """
    command = [get_bustools_binary_path(), 'inspect', path]
    p = run_executable(command, quiet=True)
    match = BUSTOOLS_INSPECT_PARSER.match(p.stdout.read())
    if not match:
        raise FileVerificationFailed(
            ('bustools inspect output could not be parsed for {}'.format(path))
        )
    if int(match.groupdict().get('count', 0)) == 0:
        raise FileVerificationFailed('{} has no BUS records'.format(path))


def validate_mtx(path):
    """Verify if the provided Matrix Market (.mtx) file is valid.

    A BUS file is considered valid when the file can be read with `scipy.io.mmread`.

    :param path: path to mtx file
    :type path: str

    :raises FileVerificationFailed: if the file failed verification
    """
    try:
        scipy.io.mmread(path)
    except ValueError:
        raise FileVerificationFailed(
            '{} is not a valid matrix market file'.format(path)
        )


VALIDATORS = {
    '.bus': validate_bus,
    '.mtx': validate_mtx,
}


def validate(path):
    """Validate a file.

    This function is a wrapper around all validation functions.
    Given a path, it chooses the correct validation function.
    This function assumes the file exists.

    :param path: path to file
    :type path: str

    :raises FileVerificationFailed: if the file failed verification
    """
    # Validation is turned off.
    if not is_validate():
        return

    ext = os.path.splitext(path)[1]
    if ext in VALIDATORS:
        VALIDATORS[ext](path)
        logger.debug('{} passed validation'.format(path))


def validate_files(pre=True, post=True):
    """Function decorator to validate input/output files.

    This function does not validate when the current run is a dry run.

    :param pre: whether to validate input files, defaults to `True`
    :type pre: bool
    :param post: whether to validate output files, defaults to `True`
    :type post: bool

    :return: wrapped function
    :rtype: function
    """

    def wrapper(func):

        @functools.wraps(func)
        def inner(*args, **kwargs):
            if not is_dry() and pre:
                for arg in list(args) + list(kwargs.values()):
                    if isinstance(arg, str) and os.path.exists(arg):
                        validate(arg)

            results = func(*args, **kwargs)

            if not is_dry() and post:
                to_check = []
                if isinstance(results, str):
                    to_check.append(results)
                if isinstance(results, (list, tuple)):
                    to_check += results
                if isinstance(results, dict):
                    to_check += list(results.values())
                for arg in to_check:
                    if isinstance(arg, str) and os.path.exists(arg):
                        validate(arg)

            return results

        return inner

    return wrapper

import functools

from ..config import is_dry


def dryable(dry_func):
    """Function decorator to set a function as dryable.

    When this decorator is applied, the provided `dry_func` will be called
    instead of the actual function when the current run is a dry run.

    :param dry_func: function to call when it is a dry run
    :type dry_func: function

    :return: wrapped function
    :rtype: function
    """

    def wrapper(func):

        @functools.wraps(func)
        def inner(*args, **kwargs):
            if not is_dry():
                return func(*args, **kwargs)
            else:
                return dry_func(*args, **kwargs)

        return inner

    return wrapper


def dummy_function(*args, **kwargs):
    """A dummy function that doesn't do anything and just returns.
    Used for making functions dryable.
    """
    return

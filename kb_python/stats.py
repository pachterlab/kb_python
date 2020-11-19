import datetime as dt
import json
import sys

from . import __version__
from .config import is_dry
from .dry import dummy_function
from .dry import dryable


class Stats:
    """Class used to collect kb run statistics.
    """

    def __init__(self):
        self.kallisto_version = None
        self.bustools_version = None
        self.start_time = None
        self.call = None
        self.commands = []
        self.runtimes = []
        self.end_time = None
        self.elapsed = None
        self.version = __version__

    def start(self):
        """Start collecting statistics.

        Sets start time, the command line call,
        and the commands array to an empty list.
        """
        self.start_time = dt.datetime.now()
        self.call = ' '.join(sys.argv)
        self.commands = []

    def command(self, command, runtime=None):
        """Report a shell command was run.

        :param command: a shell command, represented as a list
        :type command: list
        :param kwargs: additional command information
        :type kwargs: dict
        """
        cmd = ' '.join(command)
        self.commands.append(cmd)
        self.runtimes.append(runtime or 'not measured')

    def end(self):
        """End collecting statistics.
        """
        self.end_time = dt.datetime.now()
        self.elapsed = (self.end_time - self.start_time).total_seconds()

    @dryable(dummy_function)
    def save(self, path):
        """Save statistics as JSON to path.

        :param path: path to JSON
        :type path: str

        :return: path to saved JSON
        :rtype: str
        """
        if not is_dry():
            with open(path, 'w') as f:
                json.dump(self.to_dict(), f, indent=4)
        return path

    def to_dict(self):
        """Convert statistics to dictionary, so that it is easily parsed
        by the report-rendering functions.
        """
        return {
            'version': self.version,
            'start_time': self.start_time.isoformat(),
            'end_time': self.end_time.isoformat(),
            'elapsed': self.elapsed,
            'call': self.call,
            'commands': self.commands,
            'runtimes': self.runtimes,
        }


STATS = Stats()

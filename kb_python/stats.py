import datetime as dt
import sys


class Stats:
    """Class used to collect kb run statistics.
    """

    def __init__(self):
        self.kallisto_version = None
        self.bustools_version = None
        self.start_time = None
        self.call = None
        self.commands = []
        self.end_time = None
        self.elapsed = None

    def start(self):
        """Start collecting statistics.

        Sets start time, the command line call,
        and the commands array to an empty list.
        """
        self.start_time = dt.datetime.now()
        self.call = ' '.join(sys.argv)
        self.commands = []

    def command(self, command):
        """Report a shell command was run.

        :param command: a shell command, represented as a list
        :type command: list
        """
        self.commands.append(' '.join(command))

    def end(self):
        """End collecting statistics.
        """
        self.end_time = dt.datetime.now()
        self.elapsed = (self.end_time - self.start_time).total_seconds()

    def to_dict(self):
        """Convert statistics to dictionary, so that it is easily parsed
        by the report-rendering functions.
        """
        return {
            'start_time': self.start_time.isoformat(),
            'end_time': self.end_time.isoformat(),
            'elapsed': self.elapsed,
            'call': self.call,
            'commands': self.commands,
        }


STATS = Stats()

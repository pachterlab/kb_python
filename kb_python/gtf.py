import logging
import re

from .utils import open_as_text

logger = logging.getLogger(__name__)


class GTF:
    """Utility class to easily read and parse GTF files.

    :param gtf_path: path to GTF file
    :type gtf_path: str
    """
    PARSER = re.compile(
        r'''
        ^(?P<seqname>.+?)\t     # chromosome
        .*?\t                   # source
        (?P<feature>.+?)\t      # feature: transcript, exon, etc.
        (?P<start>[0-9]+?)\t    # start position (1-indexed)
        (?P<end>[0-9]+?)\t      # end position (1-indexed, inclusive)
        .*?\t                   # score
        (?P<strand>\+|-|\.)\t   # +, -, . indicating strand
        .*?\t                   # frame
        (?P<group>.*)           # groups
    ''', re.VERBOSE
    )

    GROUP_PARSER = re.compile(r'(?P<key>\S+?) "(?P<value>.+?)"')

    def __init__(self, gtf_path):
        self.gtf_path = gtf_path

    @staticmethod
    def parse_entry(line):
        """Parse a single GTF entry.

        :param line: a line in the GTF file
        :type line: str

        :return: parsed GTF information
        :rtype: dict
        """
        match = GTF.PARSER.match(line)
        if match:
            groupdict = match.groupdict()
            groupdict['start'] = int(groupdict['start'])
            groupdict['end'] = int(groupdict['end'])
            groupdict['group'] = dict(
                GTF.GROUP_PARSER.findall(groupdict.get('group', ''))
            )

            return groupdict
        return None

    def entries(self):
        """Generator that yields one GTF entry at a time.

        :return: a generator that yields a dict of the GTF entry
        :rtype: generator
        """
        with open_as_text(self.gtf_path, 'r') as f:
            for line in f:
                if line.startswith('#') or line.isspace():
                    continue

                yield GTF.parse_entry(line)

    def sort(self, out_path):
        """Sort the GTF file by chromosome, start position, line number.

        :param out_path: path to generate the sorted GTF
        :type out_path: str
        """
        to_sort = []
        with open_as_text(self.gtf_path, 'r') as f:
            position = 0
            line = f.readline()
            while line:
                if not line.startswith('#') and not line.isspace():
                    entry = GTF.parse_entry(line)
                    if entry['feature'] in ('transcript', 'exon'):
                        to_sort.append(
                            (entry['seqname'], entry['start'], position)
                        )
                position = f.tell()
                line = f.readline()
        logger.debug('Sorting {} GTF entries'.format(len(to_sort)))
        to_sort.sort()

        logger.debug('Writing sorted GTF {}'.format(out_path))
        with open_as_text(self.gtf_path, 'r') as gtf, open_as_text(out_path,
                                                                   'w') as f:
            for tup in to_sort:
                position = tup[2]
                gtf.seek(position)
                f.write(gtf.readline())

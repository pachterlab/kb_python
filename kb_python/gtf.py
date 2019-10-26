import logging
import re

logger = logging.getLogger(__name__)


class GTF:
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
        with open(self.gtf_path, 'r') as f:
            for line in f:
                if line.startswith('#') or line.isspace():
                    continue

                yield GTF.parse_entry(line)

    def sort(self, out_path):
        to_sort = []
        with open(self.gtf_path, 'r') as f:
            position = 0
            line = f.readline()
            while line:
                if not line.startswith('#') and not line.isspace():
                    entry = GTF.parse_entry(line)
                    if entry['feature'] in ('transcript', 'exon',
                                            'five_prime_utr',
                                            'three_prime_utr'):
                        to_sort.append(
                            (entry['seqname'], entry['start'], position)
                        )
                position = f.tell()
                line = f.readline()
        logger.debug('Sorting {} GTF entries'.format(len(to_sort)))
        to_sort.sort()

        logger.debug('Writing sorted GTF {}'.format(out_path))
        with open(self.gtf_path, 'r') as gtf, open(out_path, 'w') as f:
            for tup in to_sort:
                position = tup[2]
                gtf.seek(position)
                f.write(gtf.readline())

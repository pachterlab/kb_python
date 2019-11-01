import logging
import re

from .gtf import GTF
from .utils import open_as_text

logger = logging.getLogger(__name__)


class FASTA:
    """Utility class to easily read and manipulate FASTA files.

    :param fasta_path: path to FASTA file
    :type fasta_path: str
    """
    PARSER = re.compile(r'^>(?P<sequence_id>\S+)(?P<group>.*)')
    GROUP_PARSER = re.compile(r'(?P<key>\S+?):(?P<value>\S+)')
    BASEPAIRS = {
        'a': 'T',
        'A': 'T',
        'c': 'G',
        'C': 'G',
        'g': 'C',
        'G': 'C',
        't': 'A',
        'T': 'A',
        'n': 'N',
        'N': 'N',
    }

    def __init__(self, fasta_path):
        self.fasta_path = fasta_path

    @staticmethod
    def make_header(seq_id, attributes):
        """Create a correctly-formatted FASTA header with the given sequence ID
        and attributes.

        :param seq_id: sequence ID
        :type seq_id: str
        :param attributes: list of key-value pairs corresponding to attributes
                           of this sequence
        :type attributes: list

        :return: FASTA header
        :rtype: str
        """
        return '>{} {}'.format(
            seq_id, ' '.join('{}:{}'.format(k, v) for k, v in attributes)
        )

    @staticmethod
    def parse_header(line):
        """Parse information from a FASTA header.

        :param line: FASTA header line
        :type line: str

        :return: parsed information
        :rtype: dict
        """
        match = FASTA.PARSER.match(line)
        if match:
            groupdict = match.groupdict()
            groupdict['group'] = dict(
                FASTA.GROUP_PARSER.findall(groupdict.get('group', ''))
            )
            return groupdict
        return None

    @staticmethod
    def reverse_complement(sequence):
        """Get the reverse complement of the given DNA sequence.

        :param sequence: DNA sequence
        :type sequence: str

        :return: reverse complement
        :rtype: str
        """
        return ''.join(FASTA.BASEPAIRS[b] for b in reversed(sequence))

    def entries(self):
        """Generator that yields one FASTA entry (sequence ID + sequence) at a time.

        :return: a generator that yields a tuple of the FASTA entry
        :rtype: generator
        """
        with open_as_text(self.fasta_path, 'r') as f:
            info = None
            sequence = ''
            for line in f:
                if line.startswith('>'):
                    if info:
                        yield info, sequence
                        sequence = ''

                    info = FASTA.parse_header(line)
                else:
                    sequence += line.strip()

            if info:
                yield info, sequence

    def sort(self, out_path):
        """Sort the FASTA file by sequence ID.

        :param out_path: path to generate the sorted FASTA
        :type out_path: str
        """
        to_sort = []
        with open_as_text(self.fasta_path, 'r') as f:
            position = 0
            line = f.readline()

            while line:
                if line.startswith('>'):
                    to_sort.append([
                        FASTA.parse_header(line)['sequence_id'], position, None
                    ])

                position = f.tell()
                line = f.readline()

        logger.debug('Sorting {} FASTA entries'.format(len(to_sort)))
        for i in range(len(to_sort) - 1):
            to_sort[i][2] = to_sort[i + 1][1]
        to_sort.sort()

        logger.debug('Writing sorted FASTA {}'.format(out_path))
        with open_as_text(self.fasta_path,
                          'r') as fasta, open_as_text(out_path, 'w') as f:
            for tup in to_sort:
                start_position = tup[1]
                end_position = tup[2]
                fasta.seek(start_position)
                if end_position is not None:
                    while fasta.tell() < end_position:
                        f.write(fasta.readline())
                else:
                    for line in fasta:
                        f.write(line)


def generate_cdna_fasta(fasta_path, gtf_path, out_path):
    """Generate a cDNA FASTA using the genome and GTF.

    This function assumes the order in which the chromosomes appear in the
    genome FASTA is identical to the order in which they appear in the GTF.
    Additionally, the GTF must be sorted by start position.

    :param fasta_path: path to genomic FASTA file
    :type fasta_path: str
    :param gtf_path: path to GTF file
    :type gtf_path: str
    :param out_path: path to cDNA FASTA to generate
    :type out_path: str

    :return: path to generated cDNA FASTA
    :rtype: str
    """
    fasta = FASTA(fasta_path)
    gtf = GTF(gtf_path)
    gtf_entries = gtf.entries()

    with open_as_text(out_path, 'w') as f:
        previous_gtf_entry = None
        for info, sequence in fasta.entries():
            sequence_id = info['sequence_id']
            logger.debug(
                'Generating cDNA from chromosome {}'.format(sequence_id)
            )

            transcript_sequences = {}
            transcript_infos = {}
            while True:
                try:
                    gtf_entry = previous_gtf_entry if previous_gtf_entry else next(
                        gtf_entries
                    )
                except StopIteration:
                    break
                previous_gtf_entry = None
                chromosome = gtf_entry['seqname']

                if sequence_id != chromosome:
                    previous_gtf_entry = gtf_entry
                    break

                start = gtf_entry['start']
                end = gtf_entry['end']
                strand = gtf_entry['strand']
                if gtf_entry['feature'] == 'exon':
                    transcript_id = gtf_entry['group']['transcript_id']
                    transcript_version = gtf_entry['group'].get(
                        'transcript_version', None
                    )

                    transcript = '{}.{}'.format(
                        transcript_id, transcript_version
                    ) if transcript_version else transcript_id
                    if transcript not in transcript_sequences:
                        transcript_sequences[transcript] = ''
                    transcript_sequences[transcript] += sequence[start - 1:end]
                elif gtf_entry['feature'] == 'transcript':
                    transcript_id = gtf_entry['group']['transcript_id']
                    transcript_version = gtf_entry['group'].get(
                        'transcript_version', None
                    )
                    transcript = '{}.{}'.format(
                        transcript_id, transcript_version
                    ) if transcript_version else transcript_id

                    gene_id = gtf_entry['group']['gene_id']
                    gene_version = gtf_entry['group'].get('gene_version', None)
                    gene = '{}.{}'.format(
                        gene_id, gene_version
                    ) if gene_version else gene_id
                    gene_name = gtf_entry['group'].get('gene_name', '')

                    if transcript not in transcript_infos:
                        attributes = [
                            ('gene_id', gene),
                            ('gene_name', gene_name),
                            ('chr', chromosome),
                            ('start', start),
                            ('end', end),
                            ('strand', strand),
                        ]
                        transcript_infos[transcript] = attributes

            logger.debug(
                'Writing {} cDNA transcripts'.format(len(transcript_sequences))
            )
            for transcript in sorted(transcript_sequences.keys()):
                exon = transcript_sequences[transcript]
                attributes = transcript_infos[transcript]
                f.write(
                    '{}\n'.format(FASTA.make_header(transcript, attributes))
                )
                f.write(
                    '{}\n'.format(
                        exon if dict(attributes)['strand'] ==
                        '+' else FASTA.reverse_complement(exon)
                    )
                )

    return out_path


def generate_intron_fasta(fasta_path, gtf_path, out_path, flank=30):
    """Generate an intron FASTA using the genome and GTF.

    This function assumes the order in which the chromosomes appear in the
    genome FASTA is identical to the order in which they appear in the GTF.
    Additionally, the GTF must be sorted by start position.
    The intron for a specific transcript is the collection of the following:
    1. transcript - exons
    2. 5' UTR
    3. 3' UTR
    Additionally, append 30-bp (k - 1 where k = 31) flanks to each intron,
    combining sections that overlap into a single FASTA entry.

    :param fasta_path: path to genomic FASTA file
    :type fasta_path: str
    :param gtf_path: path to GTF file
    :type gtf_path: str
    :param out_path: path to intron FASTA to generate
    :type out_path: str
    :param flank: the size of intron flanks, in bases, defaults to `30`
    :type flank: int, optional

    :return: path to generated intron FASTA
    :rtype: str
    """
    fasta = FASTA(fasta_path)
    gtf = GTF(gtf_path)
    gtf_entries = gtf.entries()

    with open_as_text(out_path, 'w') as f:
        previous_gtf_entry = None
        for info, sequence in fasta.entries():
            sequence_id = info['sequence_id']
            logger.debug(
                'Generating introns from chromosome {}'.format(sequence_id)
            )

            transcript_exons = {}
            transcript_infos = {}
            while True:
                try:
                    gtf_entry = previous_gtf_entry if previous_gtf_entry else next(
                        gtf_entries
                    )
                except StopIteration:
                    break
                previous_gtf_entry = None
                chromosome = gtf_entry['seqname']

                if sequence_id != chromosome:
                    previous_gtf_entry = gtf_entry
                    break

                start = gtf_entry['start']
                end = gtf_entry['end']
                strand = gtf_entry['strand']
                if gtf_entry['feature'] == 'exon':
                    transcript_id = gtf_entry['group']['transcript_id']
                    transcript_version = gtf_entry['group'].get(
                        'transcript_version', None
                    )
                    transcript = '{}.{}'.format(
                        transcript_id, transcript_version
                    ) if transcript_version else transcript_id

                    transcript_exons.setdefault(transcript,
                                                []).append((start, end))
                elif gtf_entry['feature'] == 'transcript':
                    transcript_id = gtf_entry['group']['transcript_id']
                    transcript_version = gtf_entry['group'].get(
                        'transcript_version', None
                    )
                    transcript = '{}.{}'.format(
                        transcript_id, transcript_version
                    ) if transcript_version else transcript_id

                    gene_id = gtf_entry['group']['gene_id']
                    gene_version = gtf_entry['group'].get('gene_version', None)
                    gene = '{}.{}'.format(
                        gene_id, gene_version
                    ) if gene_version else gene_id
                    gene_name = gtf_entry['group'].get('gene_name', '')

                    if transcript not in transcript_infos:
                        attributes = [
                            ('gene_id', gene),
                            ('gene_name', gene_name),
                            ('chr', chromosome),
                            ('start', start),
                            ('end', end),
                            ('strand', strand),
                        ]
                        transcript_infos[transcript] = attributes

            for transcript in sorted(transcript_exons.keys()):
                attributes = transcript_infos[transcript]

                # Find transcript interval - all exon intervals
                attributes_dict = dict(attributes)
                transcript_interval = (
                    attributes_dict['start'], attributes_dict['end']
                )
                introns = []
                exons = list(sorted(transcript_exons[transcript]))
                if exons:
                    if exons[0][0] > transcript_interval[0]:
                        introns.append(
                            (transcript_interval[0], exons[0][0] - 1)
                        )

                    for i in range(len(exons) - 1):
                        start = exons[i][1]
                        end = exons[i + 1][0]
                        introns.append((start + 1, end - 1))

                    if exons[-1][1] < transcript_interval[1]:
                        introns.append(
                            (exons[-1][1] + 1, transcript_interval[1])
                        )
                else:
                    introns.append(transcript_interval)

                index = 1
                flank_start = None
                flank_end = None
                for start, end in introns:
                    if flank_start is None:
                        flank_start = max(start - flank, transcript_interval[0])
                    if flank_end is None or start - flank <= flank_end:
                        flank_end = min(end + flank, transcript_interval[1])
                    else:
                        intron = sequence[flank_start - 1:flank_end]
                        f.write(
                            '{}\n'.format(
                                FASTA.make_header(
                                    '{}-I.{}'.format(transcript, index),
                                    attributes
                                )
                            )
                        )
                        f.write(
                            '{}\n'.format(
                                intron if dict(attributes)['strand'] ==
                                '+' else FASTA.reverse_complement(intron)
                            )
                        )
                        index += 1
                        flank_start = max(start - flank, transcript_interval[0])
                        flank_end = min(end + flank, transcript_interval[1])
                if flank_start is not None and flank_end is not None:
                    intron = sequence[flank_start - 1:flank_end]
                    f.write(
                        '{}\n'.format(
                            FASTA.make_header(
                                '{}-I.{}'.format(transcript, index), attributes
                            )
                        )
                    )
                    f.write(
                        '{}\n'.format(
                            intron if dict(attributes)['strand'] ==
                            '+' else FASTA.reverse_complement(intron)
                        )
                    )
                    index += 1

    return out_path


def generate_spliced_fasta(fasta_path, gtf_path, out_path):
    """Generate a spliced FASTA using the genome and GTF.

    This function assumes the order in which the chromosomes appear in the
    genome FASTA is identical to the order in which they appear in the GTF.
    Additionally, the GTF must be sorted by start position.
    The spliced FASTA contains entries of length 2 * (k - 1) for k = 31,
    centered around exon-exon splice junctions (any overlapping regions are
    collapsed).

    :param fasta_path: path to genomic FASTA file
    :type fasta_path: str
    :param gtf_path: path to GTF file
    :type gtf_path: str
    :param out_path: path to spliced FASTA to generate
    :type out_path: str

    :return: path to generated spliced FASTA
    :rtype: str
    """
    pass


def generate_unspliced_fasta(fasta_path, gtf_path, out_path):
    """Generate a unspliced FASTA using the genome and GTF.

    This function assumes the order in which the chromosomes appear in the
    genome FASTA is identical to the order in which they appear in the GTF.
    Additionally, the GTF must be sorted by start position.
    The spliced FASTA contains entries of length 2 * (k - 1) for k = 31,
    centered around exon-intron splice junctions + full introns
    (any overlapping regions are collapsed).

    :param fasta_path: path to genomic FASTA file
    :type fasta_path: str
    :param gtf_path: path to GTF file
    :type gtf_path: str
    :param out_path: path to unspliced FASTA to generate
    :type out_path: str

    :return: path to generated unspliced FASTA
    :rtype: str
    """
    pass

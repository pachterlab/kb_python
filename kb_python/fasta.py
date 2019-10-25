class FASTA:
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
        return '>{} {}'.format(
            seq_id, ' '.join('{}:{}'.format(k, v) for k, v in attributes)
        )

    @staticmethod
    def parse_header(line):
        return line.strip().split(' ', 1)[0][1:]

    @staticmethod
    def reverse_complement(sequence):
        return ''.join(FASTA.BASEPAIRS[b] for b in reversed(sequence))

    def entries(self):
        with open(self.fasta_path, 'r') as f:
            sequence_id = None
            sequence = ''
            for line in f:
                if line.startswith('>'):
                    if sequence_id:
                        yield sequence_id, sequence
                        sequence = ''

                    sequence_id = FASTA.parse_header(line)
                else:
                    sequence += line.strip()

            if sequence_id:
                yield sequence_id, sequence

    def sort(self, out_path):
        to_sort = []
        with open(self.fasta_path, 'r') as f:
            position = 0
            line = f.readline()

            while line:
                if line.startswith('>'):
                    to_sort.append([FASTA.parse_header(line), position, None])

                position = f.tell()
                line = f.readline()

        for i in range(len(to_sort) - 1):
            to_sort[i][2] = to_sort[i + 1][1]
        to_sort.sort()

        with open(self.fasta_path, 'r') as fasta, open(out_path, 'w') as f:
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

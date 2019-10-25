import os
import tempfile
import uuid
from unittest import TestCase

from kb_python.fasta import FASTA
from tests.mixins import TestMixin


class TestFASTA(TestMixin, TestCase):

    def test_make_header(self):
        seq_id = 'TRANSCRIPT_ID'
        attributes = [
            ('gene_id', 'GENE_ID'),
            ('gene_name', 'GENE_NAME'),
        ]
        self.assertEqual(
            '>TRANSCRIPT_ID gene_id:GENE_ID gene_name:GENE_NAME',
            FASTA.make_header(seq_id, attributes)
        )

    def test_parse_header(self):
        header = '>transcript_id TEST'
        self.assertEqual('transcript_id', FASTA.parse_header(header))

    def test_reverse_complement(self):
        sequence = 'ATCG'
        self.assertEqual('CGAT', FASTA.reverse_complement(sequence))

    def test_entries(self):
        fasta = FASTA(self.unsorted_fasta_path)
        self.assertEqual(2, len(list(fasta.entries())))

    def test_sort(self):
        out_path = os.path.join(
            tempfile.gettempdir(), '{}.gtf'.format(uuid.uuid4())
        )
        fasta = FASTA(self.unsorted_fasta_path)
        fasta.sort(out_path)

        with open(out_path, 'r') as f, open(self.sorted_fasta_path,
                                            'r') as sorted:
            self.assertEqual(f.read(), sorted.read())

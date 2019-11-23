import os
import tempfile
import uuid
from unittest import mock, TestCase

import kb_python.fasta as fasta
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
            fasta.FASTA.make_header(seq_id, attributes)
        )

    def test_parse_header(self):
        header = '>transcript_id TEST:testing'
        self.assertEqual({
            'sequence_id': 'transcript_id',
            'group': {
                'TEST': 'testing'
            }
        }, fasta.FASTA.parse_header(header))

    def test_reverse_complement(self):
        sequence = 'ATCG'
        self.assertEqual('CGAT', fasta.FASTA.reverse_complement(sequence))

    def test_entries(self):
        fa = fasta.FASTA(self.unsorted_fasta_path)
        self.assertEqual(2, len(list(fa.entries())))

    def test_sort(self):
        out_path = os.path.join(
            tempfile.gettempdir(), '{}.fa'.format(uuid.uuid4())
        )
        fa = fasta.FASTA(self.unsorted_fasta_path)
        fa.sort(out_path)

        with open(out_path, 'r') as f, open(self.sorted_fasta_path,
                                            'r') as sorted:
            self.assertEqual(f.read(), sorted.read())

    def test_generate_kite_fasta(self):
        out_path = os.path.join(
            tempfile.gettempdir(), '{}.fa'.format(uuid.uuid4())
        )
        self.assertEqual(
            (out_path, 15),
            fasta.generate_kite_fasta(self.kite_feature_path, out_path)
        )
        with open(out_path, 'r') as f, open(self.kite_fasta_path, 'r') as fa:
            self.assertEqual(fa.read(), f.read())

    def test_generate_kite_fasta_different_length(self):
        with self.assertRaises(Exception):
            out_path = os.path.join(
                tempfile.gettempdir(), '{}.fa'.format(uuid.uuid4())
            )
            fasta.generate_kite_fasta(
                self.kite_different_feature_path, out_path
            )

    def test_generate_kite_fasta_duplicate(self):
        with self.assertRaises(Exception):
            out_path = os.path.join(
                tempfile.gettempdir(), '{}.fa'.format(uuid.uuid4())
            )
            fasta.generate_kite_fasta(
                self.kite_duplicate_feature_path, out_path
            )

    def test_generate_kite_fasta_collision(self):
        with mock.patch('kb_python.fasta.logger.warning') as warning:
            out_path = os.path.join(
                tempfile.gettempdir(), '{}.fa'.format(uuid.uuid4())
            )
            self.assertEqual((out_path, 15),
                             fasta.generate_kite_fasta(
                                 self.kite_collision_feature_path, out_path
                             ))
            warning.assert_called_once()
            with open(out_path, 'r') as f, open(self.kite_collision_fasta_path,
                                                'r') as fa:
                self.assertEqual(fa.read(), f.read())

    def test_generate_cdna_fasta(self):
        out_path = os.path.join(
            tempfile.gettempdir(), '{}.fa'.format(uuid.uuid4())
        )
        self.assertEqual(
            out_path,
            fasta.generate_cdna_fasta(
                self.sorted_fasta_path, self.sorted_gtf_path, out_path
            )
        )
        with open(out_path, 'r') as f, open(self.split_cdna_fasta_path,
                                            'r') as split:
            self.assertEqual(f.read(), split.read())

    def test_generate_intron_fasta(self):
        out_path = os.path.join(
            tempfile.gettempdir(), '{}.fa'.format(uuid.uuid4())
        )
        self.assertEqual(
            out_path,
            fasta.generate_intron_fasta(
                self.sorted_fasta_path, self.sorted_gtf_path, out_path, flank=1
            )
        )
        with open(out_path, 'r') as f, open(self.split_intron_fasta_path,
                                            'r') as split:
            self.assertEqual(f.read(), split.read())

    def test_generate_spliced_fasta(self):
        out_path = os.path.join(
            tempfile.gettempdir(), '{}.fa'.format(uuid.uuid4())
        )
        fasta.generate_spliced_fasta(
            self.sorted_fasta_path, self.sorted_gtf_path, out_path
        )

    def test_generate_unspliced_fasta(self):
        out_path = os.path.join(
            tempfile.gettempdir(), '{}.fa'.format(uuid.uuid4())
        )
        fasta.generate_unspliced_fasta(
            self.sorted_fasta_path, self.sorted_gtf_path, out_path
        )

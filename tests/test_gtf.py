import os
import uuid
from unittest import TestCase

from kb_python.gtf import GTF
from tests.mixins import TestMixin


class TestGTF(TestMixin, TestCase):

    def test_parser(self):
        line = '2\thavana\tgene\t2\t3\t.\t+\t.\tgene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC";'  # noqa
        match = GTF.PARSER.match(line)
        self.assertIsNotNone(match)
        self.assertEqual(
            {
                'seqname':
                    '2',
                'feature':
                    'gene',
                'start':
                    '2',
                'end':
                    '3',
                'strand':
                    '+',
                'group':
                    'gene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC";'  # noqa
            },
            match.groupdict()
        )

    def test_group_parser(self):
        group = 'gene_id "ENSMUSG00000102693";'
        match = GTF.GROUP_PARSER.match(group)
        self.assertIsNotNone(match)
        self.assertEqual({
            'key': 'gene_id',
            'value': 'ENSMUSG00000102693'
        }, match.groupdict())

    def test_parse_entry(self):
        line = '2\thavana\tgene\t2\t3\t.\t+\t.\tgene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC";'  # noqa
        self.assertEqual({
            'seqname': '2',
            'feature': 'gene',
            'start': 2,
            'end': 3,
            'strand': '+',
            'group': {
                'gene_id': 'ENSMUSG00000102693',
                'gene_version': '1',
                'gene_name': '4933401J01Rik',
                'gene_source': 'havana',
                'gene_biotype': 'TEC'
            }
        }, GTF.parse_entry(line))

    def test_parse_entry_with_space(self):
        line = '2\thavana\tgene\t2\t3\t.\t+\t.\tgene_id "ENSMUSG00000102693 [A]"; gene_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC";'  # noqa
        self.assertEqual({
            'seqname': '2',
            'feature': 'gene',
            'start': 2,
            'end': 3,
            'strand': '+',
            'group': {
                'gene_id': 'ENSMUSG00000102693[A]',
                'gene_version': '1',
                'gene_name': '4933401J01Rik',
                'gene_source': 'havana',
                'gene_biotype': 'TEC'
            }
        }, GTF.parse_entry(line))

    def test_entries(self):
        gtf = GTF(self.unsorted_gtf_path)
        self.assertEqual(8, len(list(gtf.entries())))

    def test_sort(self):
        out_path = os.path.join(self.temp_dir, '{}.gtf'.format(uuid.uuid4()))
        gtf = GTF(self.unsorted_gtf_path)
        self.assertEqual((out_path, {'1', '2'}), gtf.sort(out_path))

        with open(out_path, 'r') as f, open(self.sorted_gtf_path,
                                            'r') as sorted:
            self.assertEqual(f.read(), sorted.read())

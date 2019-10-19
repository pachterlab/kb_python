import os
from unittest import TestCase


class TestMixin(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.technology = '10xv2'
        cls.base_dir = os.path.dirname(os.path.abspath(__file__))
        cls.fixtures_dir = os.path.join(cls.base_dir, 'fixtures')
        cls.small_gtf_path = os.path.join(cls.fixtures_dir, 'small.gtf')
        cls.gtf_path = os.path.join(cls.fixtures_dir, 'mouse_truncated.gtf')
        cls.fasta_path = os.path.join(cls.fixtures_dir, 'mouse_truncated.fasta')
        cls.t2g_path = os.path.join(
            cls.fixtures_dir, 'transcripts_to_genes.txt'
        )
        cls.index_path = os.path.join(cls.fixtures_dir, 'mouse_truncated.idx')
        cls.whitelist_path = os.path.join(cls.fixtures_dir, 'whitelist.txt')
        cls.fastqs = [
            os.path.join(cls.fixtures_dir, 'R1.fastq'),
            os.path.join(cls.fixtures_dir, 'R2.fastq'),
        ]
        cls.txnames_path = os.path.join(cls.fixtures_dir, 'transcripts.txt')
        cls.ecmap_path = os.path.join(cls.fixtures_dir, 'matrix.ec')
        cls.bus_path = os.path.join(cls.fixtures_dir, 'output.bus')
        cls.bus_s_path = os.path.join(cls.fixtures_dir, 'output.s.bus')
        cls.bus_sc_path = os.path.join(cls.fixtures_dir, 'output.s.c.bus')
        cls.bus_scs_path = os.path.join(cls.fixtures_dir, 'output.s.c.s.bus')

        cls.counts_dir = os.path.join(cls.fixtures_dir, 'counts')
        cls.matrix_path = os.path.join(cls.counts_dir, 'genes.mtx')
        cls.barcodes_path = os.path.join(cls.counts_dir, 'genes.barcodes.txt')
        cls.genes_path = os.path.join(cls.counts_dir, 'genes.genes.txt')
        cls.loom_path = os.path.join(cls.counts_dir, 'genes.loom')
        cls.h5ad_path = os.path.join(cls.counts_dir, 'genes.h5ad')

        cls.velocity_dir = os.path.join(cls.fixtures_dir, 'velocity')
        cls.cdna_small_path = os.path.join(
            cls.velocity_dir, 'human_cdna_small.fa'
        )
        cls.cdna_small_gzip_path = os.path.join(
            cls.velocity_dir, 'human_cdna_small.fa.gz'
        )
        cls.cdna_path = os.path.join(
            cls.velocity_dir, 'human_cdna_truncated.fa'
        )
        cls.intron_path = os.path.join(
            cls.velocity_dir, 'human_intron_truncated.fa.gz'
        )

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
        cls.velocity_fastqs = [
            os.path.join(cls.velocity_dir, 'R1.fastq'),
            os.path.join(cls.velocity_dir, 'R2.fastq'),
        ]
        cls.velocity_t2g_path = os.path.join(
            cls.velocity_dir, 'transcripts_to_genes.txt'
        )
        cls.velocity_cdna_t2c_path = os.path.join(
            cls.velocity_dir, 'cdna_transcripts_to_capture.txt'
        )
        cls.velocity_intron_t2c_path = os.path.join(
            cls.velocity_dir, 'intron_transcripts_to_capture.txt'
        )
        cls.velocity_bus_path = os.path.join(cls.velocity_dir, 'output.bus')
        cls.velocity_ecmap_path = os.path.join(cls.velocity_dir, 'matrix.ec')
        cls.velocity_txnames_path = os.path.join(
            cls.velocity_dir, 'transcripts.txt'
        )
        cls.velocity_bus_scs_path = os.path.join(
            cls.velocity_dir, 'output.s.c.s.bus'
        )

        cls.gtf_dir = os.path.join(cls.fixtures_dir, 'gtf')
        cls.unsorted_gtf_path = os.path.join(cls.gtf_dir, 'not_sorted.gtf')
        cls.sorted_gtf_path = os.path.join(cls.gtf_dir, 'sorted.gtf')
        cls.gtf_t2g_path = os.path.join(cls.gtf_dir, 't2g.txt')
        cls.gtf_t2g_intron_path = os.path.join(cls.gtf_dir, 't2g_intron.txt')

        cls.fasta_dir = os.path.join(cls.fixtures_dir, 'fasta')
        cls.unsorted_fasta_path = os.path.join(cls.fasta_dir, 'not_sorted.fa')
        cls.sorted_fasta_path = os.path.join(cls.fasta_dir, 'sorted.fa')
        cls.fasta_t2c_path = os.path.join(cls.fasta_dir, 't2c.txt')
        cls.split_cdna_fasta_path = os.path.join(cls.fasta_dir, 'cdna_split.fa')
        cls.split_intron_fasta_path = os.path.join(
            cls.fasta_dir, 'intron_split.fa'
        )

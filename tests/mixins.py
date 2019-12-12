import os
from unittest import TestCase


class TestMixin(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.technology = '10xv2'
        cls.base_dir = os.path.dirname(os.path.abspath(__file__))
        cls.fixtures_dir = os.path.join(cls.base_dir, 'fixtures')
        cls.test_gz_path = os.path.join(cls.fixtures_dir, 'test.gz')
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

        cls.lamanno_dir = os.path.join(cls.fixtures_dir, 'lamanno')
        cls.cdna_path = os.path.join(cls.lamanno_dir, 'cdna.fa')
        cls.intron_path = os.path.join(cls.lamanno_dir, 'intron.fz')
        cls.lamanno_fastqs = [
            os.path.join(cls.lamanno_dir, 'R1.fastq.gz'),
            os.path.join(cls.lamanno_dir, 'R2.fastq.gz'),
        ]
        cls.lamanno_t2g_path = os.path.join(cls.lamanno_dir, 't2g.txt')
        cls.lamanno_cdna_t2c_path = os.path.join(
            cls.lamanno_dir, 'cdna_t2c.txt'
        )
        cls.lamanno_intron_t2c_path = os.path.join(
            cls.lamanno_dir, 'intron_t2c.txt'
        )
        cls.lamanno_bus_path = os.path.join(cls.lamanno_dir, 'output.bus')
        cls.lamanno_ecmap_path = os.path.join(cls.lamanno_dir, 'matrix.ec')
        cls.lamanno_txnames_path = os.path.join(
            cls.lamanno_dir, 'transcripts.txt'
        )
        cls.lamanno_bus_scs_path = os.path.join(
            cls.lamanno_dir, 'output.s.c.s.bus'
        )
        cls.lamanno_counts_dir = os.path.join(cls.lamanno_dir, 'counts')
        cls.spliced_matrix_path = os.path.join(
            cls.lamanno_counts_dir, 'spliced.mtx'
        )
        cls.spliced_barcodes_path = os.path.join(
            cls.lamanno_counts_dir, 'spliced.barcodes.txt'
        )
        cls.spliced_genes_path = os.path.join(
            cls.lamanno_counts_dir, 'spliced.genes.txt'
        )
        cls.unspliced_matrix_path = os.path.join(
            cls.lamanno_counts_dir, 'unspliced.mtx'
        )
        cls.unspliced_barcodes_path = os.path.join(
            cls.lamanno_counts_dir, 'unspliced.barcodes.txt'
        )
        cls.unspliced_genes_path = os.path.join(
            cls.lamanno_counts_dir, 'unspliced.genes.txt'
        )

        cls.gtf_dir = os.path.join(cls.fixtures_dir, 'gtf')
        cls.unsorted_gtf_path = os.path.join(cls.gtf_dir, 'not_sorted.gtf')
        cls.sorted_gtf_path = os.path.join(cls.gtf_dir, 'sorted.gtf')
        cls.gtf_t2g_path = os.path.join(cls.gtf_dir, 't2g.txt')
        cls.fasta_t2g_intron_path = os.path.join(cls.gtf_dir, 't2g_intron.txt')
        cls.gtf_t2g_intron_path = os.path.join(
            cls.gtf_dir, 't2g_intron_gtf.txt'
        )
        cls.gtf_no_exon = os.path.join(cls.gtf_dir, 'sorted_no_exon.gtf')
        cls.gtf_no_transcript = os.path.join(
            cls.gtf_dir, 'sorted_no_transcript.gtf'
        )

        cls.fasta_dir = os.path.join(cls.fixtures_dir, 'fasta')
        cls.unsorted_fasta_path = os.path.join(cls.fasta_dir, 'not_sorted.fa')
        cls.sorted_fasta_path = os.path.join(cls.fasta_dir, 'sorted.fa')
        cls.fasta_t2c_path = os.path.join(cls.fasta_dir, 't2c.txt')
        cls.split_cdna_fasta_path = os.path.join(cls.fasta_dir, 'cdna_split.fa')
        cls.split_intron_fasta_path = os.path.join(
            cls.fasta_dir, 'intron_split.fa'
        )
        cls.partial_cdna_fasta_path = os.path.join(
            cls.fasta_dir, 'cdna_partial_split.fa'
        )
        cls.partial_intron_fasta_path = os.path.join(
            cls.fasta_dir, 'intron_partial_split.fa'
        )

        cls.tcc_dir = os.path.join(cls.fixtures_dir, 'tcc')
        cls.tcc_matrix_path = os.path.join(cls.tcc_dir, 'cells_x_tcc.mtx')
        cls.tcc_barcodes_path = os.path.join(
            cls.tcc_dir, 'cells_x_tcc.barcodes.txt'
        )
        cls.tcc_ec_path = os.path.join(cls.tcc_dir, 'cells_x_tcc.ec.txt')
        cls.tcc_txnames_path = os.path.join(cls.tcc_dir, 'transcripts.txt')

        cls.kite_dir = os.path.join(cls.fixtures_dir, 'kite')
        cls.kite_feature_path = os.path.join(cls.kite_dir, 'features.tsv')
        cls.kite_collision_feature_path = os.path.join(
            cls.kite_dir, 'features_collision.tsv'
        )
        cls.kite_duplicate_feature_path = os.path.join(
            cls.kite_dir, 'features_duplicate.tsv'
        )
        cls.kite_different_feature_path = os.path.join(
            cls.kite_dir, 'features_different.tsv'
        )
        cls.kite_order_feature_path = os.path.join(
            cls.kite_dir, 'features_order.tsv'
        )
        cls.kite_fasta_path = os.path.join(cls.kite_dir, 'fasta.fa')
        cls.kite_collision_fasta_path = os.path.join(
            cls.kite_dir, 'fasta_collision.fa'
        )
        cls.kite_different_fasta_path = os.path.join(
            cls.kite_dir, 'fasta_different.fa'
        )
        cls.kite_no_mismatches_fasta_path = os.path.join(
            cls.kite_dir, 'fasta_no_mismatches.fa'
        )
        cls.kite_t2g_path = os.path.join(cls.kite_dir, 't2g.txt')
        cls.kite_map_path = os.path.join(cls.kite_dir, 'map.txt')

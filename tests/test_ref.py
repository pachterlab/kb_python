import os
import tarfile
import uuid
from unittest import mock, TestCase
from unittest.mock import call

import kb_python.ref as ref
from kb_python.config import REFERENCES_MAPPING
from tests.mixins import TestMixin


class TestRef(TestMixin, TestCase):

    def test_generate_kite_fasta(self):
        out_path = os.path.join(self.temp_dir, '{}.fa'.format(uuid.uuid4()))
        self.assertEqual(
            (out_path, 15),
            ref.generate_kite_fasta(self.kite_feature_path, out_path)
        )
        with open(out_path, 'r') as f, open(self.kite_fasta_path, 'r') as fa:
            self.assertEqual(fa.read(), f.read())

    def test_generate_kite_fasta_no_mismatches(self):
        out_path = os.path.join(self.temp_dir, '{}.fa'.format(uuid.uuid4()))
        self.assertEqual((out_path, 15),
                         ref.generate_kite_fasta(
                             self.kite_feature_path,
                             out_path,
                             no_mismatches=True
                         ))
        with open(out_path, 'r') as f, open(self.kite_no_mismatches_fasta_path,
                                            'r') as fa:
            self.assertEqual(fa.read(), f.read())

    def test_generate_kite_fasta_different_length(self):
        with mock.patch('kb_python.ref.logger.warning') as warning:
            out_path = os.path.join(self.temp_dir, '{}.fa'.format(uuid.uuid4()))
            self.assertEqual((out_path, 14),
                             ref.generate_kite_fasta(
                                 self.kite_different_feature_path, out_path
                             ))
            warning.assert_called_once()
            with open(out_path, 'r') as f, open(self.kite_different_fasta_path,
                                                'r') as fa:
                self.assertEqual(fa.read(), f.read())

    def test_generate_kite_fasta_duplicate(self):
        with self.assertRaises(Exception):
            out_path = os.path.join(self.temp_dir, '{}.fa'.format(uuid.uuid4()))
            ref.generate_kite_fasta(self.kite_duplicate_feature_path, out_path)

    def test_generate_kite_fasta_wrong_order(self):
        with self.assertRaises(Exception):
            out_path = os.path.join(self.temp_dir, '{}.fa'.format(uuid.uuid4()))
            ref.generate_kite_fasta(self.kite_order_feature_path, out_path)

    def test_generate_kite_fasta_collision(self):
        out_path = os.path.join(self.temp_dir, '{}.fa'.format(uuid.uuid4()))
        self.assertEqual(
            (out_path, 15),
            ref.generate_kite_fasta(self.kite_collision_feature_path, out_path)
        )
        with open(out_path, 'r') as f, open(self.kite_collision_fasta_path,
                                            'r') as fa:
            self.assertEqual(fa.read(), f.read())

    def test_kallisto_index(self):
        index_path = os.path.join(self.temp_dir, '{}.idx'.format(uuid.uuid4()))
        self.assertFalse(os.path.exists(index_path))
        result = ref.kallisto_index(self.fasta_path, index_path)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_split_and_index(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index:
            temp_dir = self.temp_dir
            index_prefix = os.path.join(temp_dir, 'index')
            get_temporary_filename.side_effect = [
                os.path.join(temp_dir, 'temp1'),
                os.path.join(temp_dir, 'temp2'),
                os.path.join(temp_dir, 'temp3')
            ]
            kallisto_index.side_effect = [{
                'index': '1'
            }, {
                'index': '2'
            }, {
                'index': '3'
            }]
            self.assertEqual({'indices': ['1', '2', '3']},
                             ref.split_and_index(
                                 self.fasta_path,
                                 index_prefix,
                                 temp_dir=temp_dir,
                                 n=3,
                                 k=1
                             ))
            self.assertEqual(3, kallisto_index.call_count)

    def test_create_t2g_from_fasta(self):
        t2g_path = os.path.join(self.temp_dir, '{}.txt'.format(uuid.uuid4()))
        result = ref.create_t2g_from_fasta(
            self.split_intron_fasta_path, t2g_path
        )
        with open(result['t2g'], 'r') as f, open(self.fasta_t2g_intron_path,
                                                 'r') as t2g:
            self.assertEqual(f.read(), t2g.read())

    def test_create_t2g_from_fasta_kite(self):
        t2g_path = os.path.join(self.temp_dir, '{}.txt'.format(uuid.uuid4()))
        result = ref.create_t2g_from_fasta(self.kite_fasta_path, t2g_path)
        with open(result['t2g'], 'r') as f, open(self.kite_t2g_path,
                                                 'r') as t2g:
            self.assertEqual(f.read(), t2g.read())

    def test_create_t2c(self):
        t2c_path = os.path.join(self.temp_dir, '{}.txt'.format(uuid.uuid4()))
        result = ref.create_t2c(self.unsorted_fasta_path, t2c_path)
        with open(result['t2c'], 'r') as f, open(self.fasta_t2c_path,
                                                 'r') as t2c:
            self.assertEqual(f.read(), t2c.read())

    # def test_download_reference(self):
    #     with mock.patch('kb_python.ref.download_file') as download_file:
    #         reference = REFERENCES_MAPPING['human']
    #         files = {
    #             'i': os.path.join(self.temp_dir, 'TEST.idx'),
    #             'g': os.path.join(self.temp_dir, 'TEST.txt')
    #         }
    #         temp_dir = self.temp_dir
    # 
    #         test_index_path = os.path.join(self.temp_dir, 'transcriptome.idx')
    #         test_t2g_path = os.path.join(
    #             self.temp_dir, 'transcripts_to_genes.txt'
    #         )
    #         with open(test_index_path, 'w') as index, open(test_t2g_path,
    #                                                        'w') as t2g:
    #             index.write('INDEX')
    #             t2g.write('T2G')
    #         test_tar_path = os.path.join(
    #             self.temp_dir, '{}.tar.gz'.format(uuid.uuid4())
    #         )
    #         with tarfile.open(test_tar_path, 'w:gz') as f:
    #             f.add(
    #                 test_index_path, arcname=os.path.basename(test_index_path)
    #             )
    #             f.add(test_t2g_path, arcname=os.path.basename(test_t2g_path))
    #         download_file.return_value = test_tar_path
    #         self.assertEqual(
    #             files,
    #             ref.download_reference(reference, files, temp_dir=temp_dir)
    #         )
    #         download_file.assert_called_once_with(
    #             reference.url,
    #             os.path.join(temp_dir, os.path.basename(reference.url))
    #         )
    #         with open(files['i'], 'r') as index, open(files['g'], 'r') as t2g:
    #             self.assertEqual('INDEX', index.read())
    #             self.assertEqual('T2G', t2g.read())
    # 
    # def test_download_reference_doesnt_overwrite(self):
    #     with mock.patch('kb_python.ref.os.path.exists') as exists,\
    #         mock.patch('kb_python.ref.download_file') as download_file:
    #         exists.return_value = True
    #         reference = REFERENCES_MAPPING['human']
    #         files = {
    #             'i': os.path.join(self.temp_dir, 'TEST.idx'),
    #             'g': os.path.join(self.temp_dir, 'TEST.txt')
    #         }
    #         temp_dir = self.temp_dir
    # 
    #         test_index_path = os.path.join(self.temp_dir, 'transcriptome.idx')
    #         test_t2g_path = os.path.join(
    #             self.temp_dir, 'transcripts_to_genes.txt'
    #         )
    #         with open(test_index_path, 'w') as index, open(test_t2g_path,
    #                                                        'w') as t2g:
    #             index.write('INDEX')
    #             t2g.write('T2G')
    #         test_tar_path = os.path.join(
    #             self.temp_dir, '{}.tar.gz'.format(uuid.uuid4())
    #         )
    #         with tarfile.open(test_tar_path, 'w:gz') as f:
    #             f.add(
    #                 test_index_path, arcname=os.path.basename(test_index_path)
    #             )
    #             f.add(test_t2g_path, arcname=os.path.basename(test_t2g_path))
    #         download_file.return_value = test_tar_path
    #         self.assertEqual({},
    #                          ref.download_reference(
    #                              reference, files, temp_dir=temp_dir
    #                          ))
    #         download_file.assert_not_called()
    # 
    # def test_download_reference_less_files(self):
    #     with mock.patch('kb_python.ref.download_file') as download_file:
    #         reference = REFERENCES_MAPPING['human']
    #         files = {'i': os.path.join(self.temp_dir, 'TEST.idx')}
    #         temp_dir = self.temp_dir
    # 
    #         test_index_path = os.path.join(self.temp_dir, 'transcriptome.idx')
    #         test_t2g_path = os.path.join(
    #             self.temp_dir, 'transcripts_to_genes.txt'
    #         )
    #         with open(test_index_path, 'w') as index, open(test_t2g_path,
    #                                                        'w') as t2g:
    #             index.write('INDEX')
    #             t2g.write('T2G')
    #         test_tar_path = os.path.join(
    #             self.temp_dir, '{}.tar.gz'.format(uuid.uuid4())
    #         )
    #         with tarfile.open(test_tar_path, 'w:gz') as f:
    #             f.add(
    #                 test_index_path, arcname=os.path.basename(test_index_path)
    #             )
    #             f.add(test_t2g_path, arcname=os.path.basename(test_t2g_path))
    #         download_file.return_value = test_tar_path
    #         with self.assertRaises(Exception):
    #             ref.download_reference(reference, files, temp_dir=temp_dir)
    #         download_file.assert_not_called()

    def test_decompress_file_text(self):
        with mock.patch('kb_python.ref.decompress_gzip') as decompress_gzip:
            temp_dir = self.temp_dir
            self.assertEqual(
                'textfile.txt',
                ref.decompress_file('textfile.txt', temp_dir=temp_dir)
            )
            decompress_gzip.assert_not_called()

    def test_decompress_file_gzip(self):
        with mock.patch('kb_python.ref.decompress_gzip') as decompress_gzip:
            temp_dir = self.temp_dir
            decompress_gzip.return_value = 'textfile.txt'
            self.assertEqual(
                'textfile.txt',
                ref.decompress_file('textfile.txt.gz', temp_dir=temp_dir)
            )
            decompress_gzip.assert_called_once_with(
                'textfile.txt.gz', os.path.join(temp_dir, 'textfile.txt')
            )

    def test_get_gtf_attribute_include_func(self):
        func = ref.get_gtf_attribute_include_func([{
            'a': 'b',
            'c': 'd'
        }, {
            'e': 'f'
        }])
        entry1 = mock.MagicMock()
        entry1.attributes = {'a': 'b', 'c': 'd'}
        entry2 = mock.MagicMock()
        entry2.attributes = {'a': 'b'}
        entry3 = mock.MagicMock()
        entry3.attributes = {'e': 'f'}
        entry4 = mock.MagicMock()
        entry4.attributes = {'e': 'g'}

        self.assertTrue(func(entry1))
        self.assertFalse(func(entry2))
        self.assertTrue(func(entry3))
        self.assertFalse(func(entry4))

    def test_get_gtf_attribute_exclude_func(self):
        func = ref.get_gtf_attribute_exclude_func([{
            'a': 'b',
            'c': 'd'
        }, {
            'e': 'f'
        }])
        entry1 = mock.MagicMock()
        entry1.attributes = {'a': 'b', 'c': 'd'}
        entry2 = mock.MagicMock()
        entry2.attributes = {'a': 'b'}
        entry3 = mock.MagicMock()
        entry3.attributes = {'e': 'f'}
        entry4 = mock.MagicMock()
        entry4.attributes = {'e': 'g'}

        self.assertFalse(func(entry1))
        self.assertTrue(func(entry2))
        self.assertFalse(func(entry3))
        self.assertTrue(func(entry4))

    def test_ref(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=False),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=[]):
            temp_dir = self.temp_dir
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.return_value = 'cdna'
            split_genomic_fasta_to_cdna.return_value = cdna_fasta_path
            concatenate_files.return_value = cdna_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            self.assertEqual({
                't2g': t2g_path,
                'index': index_path,
                'cdna_fasta': cdna_fasta_path,
            },
                             ref.ref(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 index_path,
                                 t2g_path,
                                 temp_dir=temp_dir
                             ))
            genes_and_transcripts_from_gtf.assert_called_once_with(
                self.gtf_path, use_version=True, filter_func=mock.ANY
            )
            create_t2g_from_fasta.assert_called_once_with(
                cdna_fasta_path, t2g_path, aa_flag=False
            )
            split_genomic_fasta_to_cdna.assert_called_once_with(
                self.fasta_path,
                'cdna',
                gene_infos,
                transcript_infos,
            )
            concatenate_files.assert_called_once_with(
                cdna_fasta_path, out_path=cdna_fasta_path
            )
            kallisto_index.assert_called_once_with(
                cdna_fasta_path,
                index_path,
                k=31,
                threads=8,
                dlist=None,
                dlist_overhang=1,
                aa=False,
                make_unique=False,
                max_ec_size=None,
                temp_dir=temp_dir
            )
            split_and_index.assert_not_called()

    def test_ref_split(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=False),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=[]):
            temp_dir = self.temp_dir
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.return_value = 'cdna'
            split_genomic_fasta_to_cdna.return_value = cdna_fasta_path
            concatenate_files.return_value = cdna_fasta_path
            split_and_index.return_value = {'indices': ['index.0', 'index.1']}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            self.assertEqual({
                't2g': t2g_path,
                'indices': ['index.0', 'index.1'],
                'cdna_fasta': cdna_fasta_path,
            },
                             ref.ref(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 index_path,
                                 t2g_path,
                                 n=2,
                                 temp_dir=temp_dir
                             ))
            genes_and_transcripts_from_gtf.assert_called_once_with(
                self.gtf_path, use_version=True, filter_func=mock.ANY
            )
            create_t2g_from_fasta.assert_called_once_with(
                cdna_fasta_path, t2g_path, aa_flag=False
            )
            split_genomic_fasta_to_cdna.assert_called_once_with(
                self.fasta_path,
                'cdna',
                gene_infos,
                transcript_infos,
            )
            concatenate_files.assert_called_once_with(
                cdna_fasta_path, out_path=cdna_fasta_path
            )
            kallisto_index.assert_not_called()
            split_and_index.assert_called_once_with(
                cdna_fasta_path, index_path, k=31, n=2, temp_dir=temp_dir
            )

    def test_ref_override_k(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=False),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=[]):
            temp_dir = self.temp_dir
            k = 999
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.return_value = 'cdna'
            split_genomic_fasta_to_cdna.return_value = cdna_fasta_path
            concatenate_files.return_value = cdna_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            self.assertEqual({
                't2g': t2g_path,
                'index': index_path,
                'cdna_fasta': cdna_fasta_path,
            },
                             ref.ref(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 index_path,
                                 t2g_path,
                                 k=k,
                                 temp_dir=temp_dir
                             ))
            genes_and_transcripts_from_gtf.assert_called_once_with(
                self.gtf_path, use_version=True, filter_func=mock.ANY
            )
            create_t2g_from_fasta.assert_called_once_with(
                cdna_fasta_path, t2g_path, aa_flag=False
            )
            split_genomic_fasta_to_cdna.assert_called_once_with(
                self.fasta_path,
                'cdna',
                gene_infos,
                transcript_infos,
            )
            concatenate_files.assert_called_once_with(
                cdna_fasta_path, out_path=cdna_fasta_path
            )
            kallisto_index.assert_called_once_with(
                cdna_fasta_path,
                index_path,
                k=k,
                threads=8,
                dlist=None,
                dlist_overhang=1,
                aa=False,
                make_unique=False,
                max_ec_size=None,
                temp_dir=temp_dir
            )
            split_and_index.assert_not_called()

    def test_ref_exists(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=True),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=[]):
            temp_dir = self.temp_dir
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.return_value = 'cdna'
            split_genomic_fasta_to_cdna.return_value = cdna_fasta_path
            concatenate_files.return_value = cdna_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            self.assertEqual({
                'index': index_path,
                't2g': t2g_path
            },
                             ref.ref(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 index_path,
                                 t2g_path,
                                 temp_dir=temp_dir
                             ))
            genes_and_transcripts_from_gtf.assert_not_called()
            split_genomic_fasta_to_cdna.assert_not_called()
            concatenate_files.assert_not_called()
            create_t2g_from_fasta.assert_called_once_with(
                cdna_fasta_path, t2g_path, aa_flag=False
            )
            kallisto_index.assert_called_once_with(
                cdna_fasta_path,
                index_path,
                k=31,
                threads=8,
                dlist=None,
                dlist_overhang=1,
                aa=False,
                make_unique=False,
                max_ec_size=None,
                temp_dir=temp_dir
            )
            split_and_index.assert_not_called()

    def test_ref_exists2(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=True),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=['a']):
            temp_dir = self.temp_dir
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.return_value = 'cdna'
            split_genomic_fasta_to_cdna.return_value = cdna_fasta_path
            concatenate_files.return_value = cdna_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            self.assertEqual({},
                             ref.ref(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 index_path,
                                 t2g_path,
                                 temp_dir=temp_dir
                             ))
            genes_and_transcripts_from_gtf.assert_not_called()
            create_t2g_from_fasta.assert_not_called()
            split_genomic_fasta_to_cdna.assert_not_called()
            concatenate_files.assert_not_called()
            kallisto_index.assert_not_called()
            split_and_index.assert_not_called()

    def test_ref_overwrite(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=True),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=['a']):
            temp_dir = self.temp_dir
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.return_value = 'cdna'
            split_genomic_fasta_to_cdna.return_value = cdna_fasta_path
            concatenate_files.return_value = cdna_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            self.assertEqual({
                't2g': t2g_path,
                'index': index_path,
                'cdna_fasta': cdna_fasta_path,
            },
                             ref.ref(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 index_path,
                                 t2g_path,
                                 temp_dir=temp_dir,
                                 overwrite=True
                             ))
            genes_and_transcripts_from_gtf.assert_called_once_with(
                self.gtf_path, use_version=True, filter_func=mock.ANY
            )
            create_t2g_from_fasta.assert_called_once_with(
                cdna_fasta_path, t2g_path, aa_flag=False
            )
            split_genomic_fasta_to_cdna.assert_called_once_with(
                self.fasta_path,
                'cdna',
                gene_infos,
                transcript_infos,
            )
            concatenate_files.assert_called_once_with(
                cdna_fasta_path, out_path=cdna_fasta_path
            )
            kallisto_index.assert_called_once_with(
                cdna_fasta_path,
                index_path,
                k=31,
                threads=8,
                dlist=None,
                dlist_overhang=1,
                aa=False,
                make_unique=False,
                max_ec_size=None,
                temp_dir=temp_dir
            )
            split_and_index.assert_not_called()

    def test_ref_kite_odd(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.generate_kite_fasta') as generate_kite_fasta,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob') as glob:
            glob.return_value = []
            temp_dir = self.temp_dir
            feature_path = mock.MagicMock()
            fasta_path = mock.MagicMock()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            decompress_file.return_value = feature_path
            generate_kite_fasta.return_value = fasta_path, 1
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            kallisto_index.return_value = {'index': index_path}

            self.assertEqual({
                'fasta': fasta_path,
                't2g': t2g_path,
                'index': index_path,
            },
                             ref.ref_kite(
                                 feature_path,
                                 fasta_path,
                                 index_path,
                                 t2g_path,
                                 temp_dir=temp_dir
                             ))
            decompress_file.assert_called_once_with(
                feature_path, temp_dir=temp_dir
            )
            generate_kite_fasta.assert_called_once_with(
                feature_path, fasta_path, no_mismatches=False
            )
            create_t2g_from_fasta.assert_called_once_with(fasta_path, t2g_path)
            kallisto_index.assert_called_once_with(
                fasta_path, index_path, k=1, threads=8, temp_dir=temp_dir
            )
            split_and_index.assert_not_called()

    def test_ref_kite_split(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.generate_kite_fasta') as generate_kite_fasta,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob') as glob:
            glob.return_value = []
            temp_dir = self.temp_dir
            feature_path = mock.MagicMock()
            fasta_path = mock.MagicMock()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            decompress_file.return_value = feature_path
            generate_kite_fasta.return_value = fasta_path, 1
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            split_and_index.return_value = {'indices': ['index.0', 'index.1']}

            self.assertEqual({
                'fasta': fasta_path,
                't2g': t2g_path,
                'indices': ['index.0', 'index.1'],
            },
                             ref.ref_kite(
                                 feature_path,
                                 fasta_path,
                                 index_path,
                                 t2g_path,
                                 n=2,
                                 temp_dir=temp_dir
                             ))
            decompress_file.assert_called_once_with(
                feature_path, temp_dir=temp_dir
            )
            generate_kite_fasta.assert_called_once_with(
                feature_path, fasta_path, no_mismatches=False
            )
            create_t2g_from_fasta.assert_called_once_with(fasta_path, t2g_path)
            split_and_index.assert_called_once_with(
                fasta_path, index_path, k=1, n=2, temp_dir=temp_dir
            )
            kallisto_index.assert_not_called()

    def test_ref_kite_even(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.generate_kite_fasta') as generate_kite_fasta,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob') as glob:
            glob.return_value = []
            temp_dir = self.temp_dir
            feature_path = mock.MagicMock()
            fasta_path = mock.MagicMock()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            decompress_file.return_value = feature_path
            generate_kite_fasta.return_value = fasta_path, 2
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            kallisto_index.return_value = {'index': index_path}

            self.assertEqual({
                'fasta': fasta_path,
                't2g': t2g_path,
                'index': index_path,
            },
                             ref.ref_kite(
                                 feature_path,
                                 fasta_path,
                                 index_path,
                                 t2g_path,
                                 temp_dir=temp_dir
                             ))
            decompress_file.assert_called_once_with(
                feature_path, temp_dir=temp_dir
            )
            generate_kite_fasta.assert_called_once_with(
                feature_path, fasta_path, no_mismatches=False
            )
            create_t2g_from_fasta.assert_called_once_with(fasta_path, t2g_path)
            kallisto_index.assert_called_once_with(
                fasta_path, index_path, k=1, threads=8, temp_dir=temp_dir
            )

    def test_ref_kite_override_k(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.generate_kite_fasta') as generate_kite_fasta,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob') as glob:
            glob.return_value = []
            k = 999
            temp_dir = self.temp_dir
            feature_path = mock.MagicMock()
            fasta_path = mock.MagicMock()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            decompress_file.return_value = feature_path
            generate_kite_fasta.return_value = fasta_path, 2
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            kallisto_index.return_value = {'index': index_path}

            self.assertEqual({
                'fasta': fasta_path,
                't2g': t2g_path,
                'index': index_path,
            },
                             ref.ref_kite(
                                 feature_path,
                                 fasta_path,
                                 index_path,
                                 t2g_path,
                                 k=k,
                                 temp_dir=temp_dir
                             ))
            decompress_file.assert_called_once_with(
                feature_path, temp_dir=temp_dir
            )
            generate_kite_fasta.assert_called_once_with(
                feature_path, fasta_path, no_mismatches=False
            )
            create_t2g_from_fasta.assert_called_once_with(fasta_path, t2g_path)
            kallisto_index.assert_called_once_with(
                fasta_path, index_path, k=k, threads=8, temp_dir=temp_dir
            )

    def test_ref_kite_doesnt_overwrite(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.generate_kite_fasta') as generate_kite_fasta,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob') as glob:
            glob.return_value = ['index']
            temp_dir = self.temp_dir
            feature_path = mock.MagicMock()
            fasta_path = mock.MagicMock()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            decompress_file.return_value = feature_path
            generate_kite_fasta.return_value = fasta_path, 1
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}

            self.assertEqual({
                'fasta': fasta_path,
                't2g': t2g_path,
            },
                             ref.ref_kite(
                                 feature_path,
                                 fasta_path,
                                 index_path,
                                 t2g_path,
                                 temp_dir=temp_dir
                             ))
            decompress_file.assert_called_once_with(
                feature_path, temp_dir=temp_dir
            )
            generate_kite_fasta.assert_called_once_with(
                feature_path, fasta_path, no_mismatches=False
            )
            create_t2g_from_fasta.assert_called_once_with(fasta_path, t2g_path)
            kallisto_index.assert_not_called()

    def test_ref_kite_overwrite(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.generate_kite_fasta') as generate_kite_fasta,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob') as glob:
            glob.return_value = ['index']
            temp_dir = self.temp_dir
            feature_path = mock.MagicMock()
            fasta_path = mock.MagicMock()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            decompress_file.return_value = feature_path
            generate_kite_fasta.return_value = fasta_path, 1
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            kallisto_index.return_value = {'index': index_path}

            self.assertEqual({
                'fasta': fasta_path,
                't2g': t2g_path,
                'index': index_path,
            },
                             ref.ref_kite(
                                 feature_path,
                                 fasta_path,
                                 index_path,
                                 t2g_path,
                                 temp_dir=temp_dir,
                                 overwrite=True
                             ))
            decompress_file.assert_called_once_with(
                feature_path, temp_dir=temp_dir
            )
            generate_kite_fasta.assert_called_once_with(
                feature_path, fasta_path, no_mismatches=False
            )
            create_t2g_from_fasta.assert_called_once_with(fasta_path, t2g_path)
            kallisto_index.assert_called_once_with(
                fasta_path, index_path, k=1, threads=8, temp_dir=temp_dir
            )

    def test_ref_nac(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_nascent') as split_genomic_fasta_to_nascent,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=False),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=[]):
            temp_dir = self.temp_dir
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            intron_fasta_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.side_effect = [
                'cdna', 'cdna_t2c', 'intron', 'intron_t2c', 'combined'
            ]
            split_genomic_fasta_to_cdna.return_value = 'cdna'
            split_genomic_fasta_to_nascent.return_value = 'intron'
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': 'cdna_t2c'
            }, {
                't2c': 'intron_t2c'
            }]
            concatenate_files.side_effect = [
                cdna_fasta_path, cdna_t2c_path, intron_fasta_path,
                intron_t2c_path, combined_path
            ]
            self.assertEqual({
                't2g': t2g_path,
                'cdna_fasta': cdna_fasta_path,
                'intron_fasta': intron_fasta_path,
                'cdna_t2c': cdna_t2c_path,
                'intron_t2c': intron_t2c_path,
                'index': index_path,
            },
                             ref.ref_nac(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 intron_fasta_path,
                                 index_path,
                                 t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 temp_dir=temp_dir
                             ))
            genes_and_transcripts_from_gtf.assert_called_once_with(
                self.gtf_path, use_version=True, filter_func=mock.ANY
            )
            create_t2g_from_fasta.assert_called_once_with(
                combined_path, t2g_path
            )
            split_genomic_fasta_to_cdna.assert_called_once_with(
                self.fasta_path, 'cdna', gene_infos, transcript_infos
            )
            split_genomic_fasta_to_nascent.assert_called_once_with(
                self.fasta_path,
                'intron',
                gene_infos
            )
            self.assertEqual(2, create_t2c.call_count)
            create_t2c.assert_has_calls([
                call('cdna', 'cdna_t2c'),
                call('intron', 'intron_t2c')
            ])
            self.assertEqual(5, concatenate_files.call_count)
            concatenate_files.assert_has_calls([
                call('cdna', out_path=cdna_fasta_path),
                call('cdna_t2c', out_path=cdna_t2c_path),
                call('intron', out_path=intron_fasta_path),
                call('intron_t2c', out_path=intron_t2c_path),
                call(
                    cdna_fasta_path,
                    intron_fasta_path,
                    out_path='combined',
                )
            ])
            kallisto_index.assert_called_once_with(
                combined_path,
                index_path,
                k=31,
                threads=8,
                dlist=None,
                dlist_overhang=1,
                make_unique=False,
                max_ec_size=None,
                temp_dir=temp_dir
            )
            split_and_index.assert_not_called()

    def test_ref_lamanno_split_2(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_intron') as split_genomic_fasta_to_intron,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=False),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=[]):
            temp_dir = self.temp_dir
            index_path = 'index'
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            intron_fasta_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.side_effect = [
                'cdna', 'cdna_t2c', 'intron', 'intron_t2c', 'combined'
            ]
            split_genomic_fasta_to_cdna.return_value = 'cdna'
            split_genomic_fasta_to_intron.return_value = 'intron'
            kallisto_index.side_effect = [{
                'index': 'index_cdna'
            }, {
                'index': 'index_intron'
            }]
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': 'cdna_t2c'
            }, {
                't2c': 'intron_t2c'
            }]
            concatenate_files.side_effect = [
                cdna_fasta_path, cdna_t2c_path, intron_fasta_path,
                intron_t2c_path, combined_path
            ]
            self.assertEqual({
                't2g': t2g_path,
                'cdna_fasta': cdna_fasta_path,
                'intron_fasta': intron_fasta_path,
                'cdna_t2c': cdna_t2c_path,
                'intron_t2c': intron_t2c_path,
                'indices': ['index_cdna', 'index_intron'],
            },
                             ref.ref_lamanno(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 intron_fasta_path,
                                 index_path,
                                 t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 temp_dir=temp_dir,
                                 n=2
                             ))
            genes_and_transcripts_from_gtf.assert_called_once_with(
                self.gtf_path, use_version=True, filter_func=mock.ANY
            )
            create_t2g_from_fasta.assert_called_once_with(
                combined_path, t2g_path
            )
            split_genomic_fasta_to_cdna.assert_called_once_with(
                self.fasta_path, 'cdna', gene_infos, transcript_infos
            )
            split_genomic_fasta_to_intron.assert_called_once_with(
                self.fasta_path,
                'intron',
                gene_infos,
                transcript_infos,
                flank=30
            )
            self.assertEqual(2, create_t2c.call_count)
            create_t2c.assert_has_calls([
                call('cdna', 'cdna_t2c'),
                call('intron', 'intron_t2c')
            ])
            self.assertEqual(5, concatenate_files.call_count)
            concatenate_files.assert_has_calls([
                call('cdna', out_path=cdna_fasta_path),
                call('cdna_t2c', out_path=cdna_t2c_path),
                call('intron', out_path=intron_fasta_path),
                call('intron_t2c', out_path=intron_t2c_path),
                call(
                    cdna_fasta_path,
                    intron_fasta_path,
                    out_path='combined',
                )
            ])
            self.assertEqual(2, kallisto_index.call_count)
            split_and_index.assert_not_called()

    def test_ref_lamanno_split_3(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_intron') as split_genomic_fasta_to_intron,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=False),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=[]):
            temp_dir = self.temp_dir
            index_path = 'index'
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            intron_fasta_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.side_effect = [
                'cdna', 'cdna_t2c', 'intron', 'intron_t2c', 'combined'
            ]
            split_genomic_fasta_to_cdna.return_value = 'cdna'
            split_genomic_fasta_to_intron.return_value = 'intron'
            kallisto_index.return_value = {'index': 'index_cdna'}
            split_and_index.return_value = {
                'indices': ['index_intron.0', 'index_intron.1']
            }
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': 'cdna_t2c'
            }, {
                't2c': 'intron_t2c'
            }]
            concatenate_files.side_effect = [
                cdna_fasta_path, cdna_t2c_path, intron_fasta_path,
                intron_t2c_path, combined_path
            ]
            self.assertEqual({
                't2g': t2g_path,
                'cdna_fasta': cdna_fasta_path,
                'intron_fasta': intron_fasta_path,
                'cdna_t2c': cdna_t2c_path,
                'intron_t2c': intron_t2c_path,
                'indices': ['index_cdna', 'index_intron.0', 'index_intron.1'],
            },
                             ref.ref_lamanno(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 intron_fasta_path,
                                 index_path,
                                 t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 temp_dir=temp_dir,
                                 n=3
                             ))
            genes_and_transcripts_from_gtf.assert_called_once_with(
                self.gtf_path, use_version=True, filter_func=mock.ANY
            )
            create_t2g_from_fasta.assert_called_once_with(
                combined_path, t2g_path
            )
            split_genomic_fasta_to_cdna.assert_called_once_with(
                self.fasta_path, 'cdna', gene_infos, transcript_infos
            )
            split_genomic_fasta_to_intron.assert_called_once_with(
                self.fasta_path,
                'intron',
                gene_infos,
                transcript_infos,
                flank=30
            )
            self.assertEqual(2, create_t2c.call_count)
            create_t2c.assert_has_calls([
                call('cdna', 'cdna_t2c'),
                call('intron', 'intron_t2c')
            ])
            self.assertEqual(5, concatenate_files.call_count)
            concatenate_files.assert_has_calls([
                call('cdna', out_path=cdna_fasta_path),
                call('cdna_t2c', out_path=cdna_t2c_path),
                call('intron', out_path=intron_fasta_path),
                call('intron_t2c', out_path=intron_t2c_path),
                call(
                    cdna_fasta_path,
                    intron_fasta_path,
                    out_path='combined',
                )
            ])
            kallisto_index.assert_called_once_with(
                cdna_fasta_path, 'index_cdna', k=31, temp_dir=temp_dir
            )
            split_and_index.assert_called_once_with(
                intron_fasta_path, 'index_intron', n=2, k=31, temp_dir=temp_dir
            )

    def test_ref_nac_override_k(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_nascent') as split_genomic_fasta_to_nascent,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=False),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=[]):
            temp_dir = self.temp_dir
            k = 999
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            intron_fasta_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.side_effect = [
                'cdna', 'cdna_t2c', 'intron', 'intron_t2c', 'combined'
            ]
            split_genomic_fasta_to_cdna.return_value = 'cdna'
            split_genomic_fasta_to_nascent.return_value = 'intron'
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': 'cdna_t2c'
            }, {
                't2c': 'intron_t2c'
            }]
            concatenate_files.side_effect = [
                cdna_fasta_path, cdna_t2c_path, intron_fasta_path,
                intron_t2c_path, combined_path
            ]
            self.assertEqual({
                't2g': t2g_path,
                'cdna_fasta': cdna_fasta_path,
                'intron_fasta': intron_fasta_path,
                'cdna_t2c': cdna_t2c_path,
                'intron_t2c': intron_t2c_path,
                'index': index_path,
            },
                             ref.ref_nac(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 intron_fasta_path,
                                 index_path,
                                 t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 k=k,
                                 temp_dir=temp_dir
                             ))
            genes_and_transcripts_from_gtf.assert_called_once_with(
                self.gtf_path, use_version=True, filter_func=mock.ANY
            )
            create_t2g_from_fasta.assert_called_once_with(
                combined_path, t2g_path
            )
            split_genomic_fasta_to_cdna.assert_called_once_with(
                self.fasta_path, 'cdna', gene_infos, transcript_infos
            )
            split_genomic_fasta_to_nascent.assert_called_once_with(
                self.fasta_path,
                'intron',
                gene_infos
            )
            self.assertEqual(2, create_t2c.call_count)
            create_t2c.assert_has_calls([
                call('cdna', 'cdna_t2c'),
                call('intron', 'intron_t2c')
            ])
            self.assertEqual(5, concatenate_files.call_count)
            concatenate_files.assert_has_calls([
                call('cdna', out_path=cdna_fasta_path),
                call('cdna_t2c', out_path=cdna_t2c_path),
                call('intron', out_path=intron_fasta_path),
                call('intron_t2c', out_path=intron_t2c_path),
                call(
                    cdna_fasta_path,
                    intron_fasta_path,
                    out_path='combined',
                )
            ])
            kallisto_index.assert_called_once_with(
                combined_path,
                index_path,
                k=k,
                threads=8,
                dlist=None,
                dlist_overhang=1,
                make_unique=False,
                max_ec_size=None,
                temp_dir=temp_dir
            )
            split_and_index.assert_not_called()

    def test_ref_nac_exists(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_intron') as split_genomic_fasta_to_intron,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=True),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=[]):
            temp_dir = self.temp_dir
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            intron_fasta_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.return_value = 'combined'
            split_genomic_fasta_to_cdna.return_value = 'cdna'
            split_genomic_fasta_to_intron.return_value = 'intron'
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': 'cdna_t2c'
            }, {
                't2c': 'intron_t2c'
            }]
            concatenate_files.return_value = combined_path
            self.assertEqual({
                't2g': t2g_path,
                'index': index_path,
            },
                             ref.ref_nac(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 intron_fasta_path,
                                 index_path,
                                 t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 temp_dir=temp_dir
                             ))
            genes_and_transcripts_from_gtf.assert_not_called()
            create_t2g_from_fasta.assert_called_once_with(
                combined_path, t2g_path
            )
            split_genomic_fasta_to_cdna.assert_not_called()
            split_genomic_fasta_to_intron.assert_not_called()
            create_t2c.assert_not_called()
            concatenate_files.assert_called_once_with(
                cdna_fasta_path,
                intron_fasta_path,
                out_path='combined',
            )
            kallisto_index.assert_called_once_with(
                combined_path,
                index_path,
                k=31,
                threads=8,
                dlist=None,
                dlist_overhang=1,
                make_unique=False,
                max_ec_size=None,
                temp_dir=temp_dir
            )
            split_and_index.assert_not_called()

    def test_ref_nac_exists2(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_intron') as split_genomic_fasta_to_intron,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=True),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=['a']):
            temp_dir = self.temp_dir
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            intron_fasta_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.return_value = 'combined'
            split_genomic_fasta_to_cdna.return_value = 'cdna'
            split_genomic_fasta_to_intron.return_value = 'intron'
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': 'cdna_t2c'
            }, {
                't2c': 'intron_t2c'
            }]
            concatenate_files.return_value = combined_path
            self.assertEqual({},
                             ref.ref_nac(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 intron_fasta_path,
                                 index_path,
                                 t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 temp_dir=temp_dir
                             ))
            genes_and_transcripts_from_gtf.assert_not_called()
            create_t2g_from_fasta.assert_not_called()
            split_genomic_fasta_to_cdna.assert_not_called()
            split_genomic_fasta_to_intron.assert_not_called()
            create_t2c.assert_not_called()
            concatenate_files.assert_not_called()
            kallisto_index.assert_not_called()
            split_and_index.assert_not_called()

    def test_ref_nac_overwrite(self):
        with mock.patch('kb_python.ref.get_temporary_filename') as get_temporary_filename,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.ngs.gtf.genes_and_transcripts_from_gtf') as genes_and_transcripts_from_gtf,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_cdna') as split_genomic_fasta_to_cdna,\
            mock.patch('kb_python.ref.ngs.fasta.split_genomic_fasta_to_nascent') as split_genomic_fasta_to_nascent,\
            mock.patch('kb_python.ref.ngs.utils.all_exists', return_value=True),\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.split_and_index') as split_and_index,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.glob.glob', return_value=['a']):
            temp_dir = self.temp_dir
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            intron_fasta_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            gene_infos = mock.MagicMock()
            transcript_infos = mock.MagicMock()
            genes_and_transcripts_from_gtf.return_value = (
                gene_infos, transcript_infos
            )
            get_temporary_filename.side_effect = [
                'cdna', 'cdna_t2c', 'intron', 'intron_t2c', 'combined'
            ]
            split_genomic_fasta_to_cdna.return_value = 'cdna'
            split_genomic_fasta_to_nascent.return_value = 'intron'
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': 'cdna_t2c'
            }, {
                't2c': 'intron_t2c'
            }]
            concatenate_files.side_effect = [
                cdna_fasta_path, cdna_t2c_path, intron_fasta_path,
                intron_t2c_path, combined_path
            ]
            self.assertEqual({
                't2g': t2g_path,
                'cdna_fasta': cdna_fasta_path,
                'intron_fasta': intron_fasta_path,
                'cdna_t2c': cdna_t2c_path,
                'intron_t2c': intron_t2c_path,
                'index': index_path,
            },
                             ref.ref_nac(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 intron_fasta_path,
                                 index_path,
                                 t2g_path,
                                 cdna_t2c_path,
                                 intron_t2c_path,
                                 temp_dir=temp_dir,
                                 overwrite=True
                             ))
            genes_and_transcripts_from_gtf.assert_called_once_with(
                self.gtf_path, use_version=True, filter_func=mock.ANY
            )
            create_t2g_from_fasta.assert_called_once_with(
                combined_path, t2g_path
            )
            split_genomic_fasta_to_cdna.assert_called_once_with(
                self.fasta_path, 'cdna', gene_infos, transcript_infos
            )
            split_genomic_fasta_to_nascent.assert_called_once_with(
                self.fasta_path,
                'intron',
                gene_infos
            )
            self.assertEqual(2, create_t2c.call_count)
            create_t2c.assert_has_calls([
                call('cdna', 'cdna_t2c'),
                call('intron', 'intron_t2c')
            ])
            self.assertEqual(5, concatenate_files.call_count)
            concatenate_files.assert_has_calls([
                call('cdna', out_path=cdna_fasta_path),
                call('cdna_t2c', out_path=cdna_t2c_path),
                call('intron', out_path=intron_fasta_path),
                call('intron_t2c', out_path=intron_t2c_path),
                call(
                    cdna_fasta_path,
                    intron_fasta_path,
                    out_path='combined',
                )
            ])
            kallisto_index.assert_called_once_with(
                combined_path,
                index_path,
                k=31,
                threads=8,
                dlist=None,
                dlist_overhang=1,
                make_unique=False,
                max_ec_size=None,
                temp_dir=temp_dir
            )
            split_and_index.assert_not_called()

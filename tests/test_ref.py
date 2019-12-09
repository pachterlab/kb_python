import os
import tarfile
import tempfile
import uuid
from unittest import mock, TestCase
from unittest.mock import call

import kb_python.ref as ref
from kb_python.config import REFERENCES_MAPPING
from kb_python.constants import (
    SORTED_FASTA_FILENAME,
    SORTED_GTF_FILENAME,
)
from tests.mixins import TestMixin


class TestRef(TestMixin, TestCase):

    def test_sort_gtf(self):
        with mock.patch('kb_python.ref.GTF') as GTF:
            out_path = 'test'
            self.assertEqual(
                GTF().sort.return_value,
                ref.sort_gtf(self.unsorted_gtf_path, out_path)
            )
            GTF().sort.assert_called_once_with(out_path)

    def test_sort_fasta(self):
        with mock.patch('kb_python.ref.FASTA') as FASTA:
            out_path = 'test'
            self.assertEqual(
                FASTA().sort.return_value,
                ref.sort_fasta(self.unsorted_fasta_path, out_path)
            )
            FASTA().sort.assert_called_once_with(out_path)

    def test_check_chromosomes(self):
        with mock.patch('kb_python.ref.logger.warning') as warning:
            fasta_chromosomes = gtf_chromosomes = chromosomes = {'1', '2'}
            self.assertEqual(
                chromosomes,
                ref.check_chromosomes(fasta_chromosomes, gtf_chromosomes)
            )
            warning.assert_not_called()

    def test_check_chromosomes_fasta_unique(self):
        with mock.patch('kb_python.ref.logger.warning') as warning:
            fasta_chromosomes = {'1', '2'}
            gtf_chromosomes = chromosomes = {'2'}
            self.assertEqual(
                chromosomes,
                ref.check_chromosomes(fasta_chromosomes, gtf_chromosomes)
            )
            warning.assert_called_once()

    def test_check_chromosomes_gtf_unique(self):
        with mock.patch('kb_python.ref.logger.warning') as warning:
            fasta_chromosomes = chromosomes = {'2'}
            gtf_chromosomes = {'1', '2'}
            self.assertEqual(
                chromosomes,
                ref.check_chromosomes(fasta_chromosomes, gtf_chromosomes)
            )
            warning.assert_called_once()

    def test_check_chromosomes_both_unique(self):
        with mock.patch('kb_python.ref.logger.warning') as warning:
            fasta_chromosomes = {'1'}
            gtf_chromosomes = {'2'}
            chromosomes = set()
            self.assertEqual(
                chromosomes,
                ref.check_chromosomes(fasta_chromosomes, gtf_chromosomes)
            )
            self.assertEqual(2, warning.call_count)

    def test_kallisto_index(self):
        index_path = os.path.join(
            tempfile.gettempdir(), '{}.idx'.format(uuid.uuid4())
        )
        self.assertFalse(os.path.exists(index_path))
        result = ref.kallisto_index(self.fasta_path, index_path)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

    def test_create_t2g_from_fasta(self):
        t2g_path = os.path.join(
            tempfile.gettempdir(), '{}.txt'.format(uuid.uuid4())
        )
        result = ref.create_t2g_from_fasta(
            self.split_intron_fasta_path, t2g_path
        )
        with open(result['t2g'], 'r') as f, open(self.fasta_t2g_intron_path,
                                                 'r') as t2g:
            self.assertEqual(f.read(), t2g.read())

    def test_create_t2g_from_fasta_kite(self):
        t2g_path = os.path.join(
            tempfile.gettempdir(), '{}.txt'.format(uuid.uuid4())
        )
        result = ref.create_t2g_from_fasta(self.kite_fasta_path, t2g_path)
        with open(result['t2g'], 'r') as f, open(self.kite_t2g_path,
                                                 'r') as t2g:
            self.assertEqual(f.read(), t2g.read())

    def test_create_t2g_from_gtf(self):
        t2g_path = os.path.join(
            tempfile.gettempdir(), '{}.txt'.format(uuid.uuid4())
        )
        result = ref.create_t2g_from_gtf(self.unsorted_gtf_path, t2g_path)
        with open(result['t2g'], 'r') as f, open(self.gtf_t2g_path, 'r') as t2g:
            self.assertEqual(f.read(), t2g.read())

    def test_create_t2g_from_gtf_with_intron(self):
        t2g_path = os.path.join(
            tempfile.gettempdir(), '{}.txt'.format(uuid.uuid4())
        )
        result = ref.create_t2g_from_gtf(
            self.unsorted_gtf_path, t2g_path, intron=True
        )
        with open(result['t2g'], 'r') as f, open(self.gtf_t2g_intron_path,
                                                 'r') as t2g:
            self.assertEqual(f.read(), t2g.read())

    def test_create_t2c(self):
        t2c_path = os.path.join(
            tempfile.gettempdir(), '{}.txt'.format(uuid.uuid4())
        )
        result = ref.create_t2c(self.unsorted_fasta_path, t2c_path)
        with open(result['t2c'], 'r') as f, open(self.fasta_t2c_path,
                                                 'r') as t2c:
            self.assertEqual(f.read(), t2c.read())

    def test_download_reference(self):
        with mock.patch('kb_python.ref.download_file') as download_file:
            reference = REFERENCES_MAPPING['human']
            files = {
                'i': os.path.join(tempfile.mkdtemp(), 'TEST.idx'),
                'g': os.path.join(tempfile.mkdtemp(), 'TEST.txt')
            }
            temp_dir = tempfile.mkdtemp()

            test_index_path = os.path.join(
                tempfile.mkdtemp(), 'transcriptome.idx'
            )
            test_t2g_path = os.path.join(
                tempfile.mkdtemp(), 'transcripts_to_genes.txt'
            )
            with open(test_index_path, 'w') as index, open(test_t2g_path,
                                                           'w') as t2g:
                index.write('INDEX')
                t2g.write('T2G')
            test_tar_path = os.path.join(
                tempfile.gettempdir(), '{}.tar.gz'.format(uuid.uuid4())
            )
            with tarfile.open(test_tar_path, 'w:gz') as f:
                f.add(
                    test_index_path, arcname=os.path.basename(test_index_path)
                )
                f.add(test_t2g_path, arcname=os.path.basename(test_t2g_path))
            download_file.return_value = test_tar_path
            self.assertEqual(
                files,
                ref.download_reference(reference, files, temp_dir=temp_dir)
            )
            download_file.assert_called_once_with(
                reference.url,
                os.path.join(temp_dir, os.path.basename(reference.url))
            )
            with open(files['i'], 'r') as index, open(files['g'], 'r') as t2g:
                self.assertEqual('INDEX', index.read())
                self.assertEqual('T2G', t2g.read())

    def test_download_reference_doesnt_overwrite(self):
        with mock.patch('kb_python.ref.os.path.exists') as exists,\
            mock.patch('kb_python.ref.download_file') as download_file:
            exists.return_value = True
            reference = REFERENCES_MAPPING['human']
            files = {
                'i': os.path.join(tempfile.mkdtemp(), 'TEST.idx'),
                'g': os.path.join(tempfile.mkdtemp(), 'TEST.txt')
            }
            temp_dir = tempfile.mkdtemp()

            test_index_path = os.path.join(
                tempfile.mkdtemp(), 'transcriptome.idx'
            )
            test_t2g_path = os.path.join(
                tempfile.mkdtemp(), 'transcripts_to_genes.txt'
            )
            with open(test_index_path, 'w') as index, open(test_t2g_path,
                                                           'w') as t2g:
                index.write('INDEX')
                t2g.write('T2G')
            test_tar_path = os.path.join(
                tempfile.gettempdir(), '{}.tar.gz'.format(uuid.uuid4())
            )
            with tarfile.open(test_tar_path, 'w:gz') as f:
                f.add(
                    test_index_path, arcname=os.path.basename(test_index_path)
                )
                f.add(test_t2g_path, arcname=os.path.basename(test_t2g_path))
            download_file.return_value = test_tar_path
            self.assertEqual({},
                             ref.download_reference(
                                 reference, files, temp_dir=temp_dir
                             ))
            download_file.assert_not_called()

    def test_download_reference_less_files(self):
        with mock.patch('kb_python.ref.download_file') as download_file:
            reference = REFERENCES_MAPPING['human']
            files = {'i': os.path.join(tempfile.mkdtemp(), 'TEST.idx')}
            temp_dir = tempfile.mkdtemp()

            test_index_path = os.path.join(
                tempfile.mkdtemp(), 'transcriptome.idx'
            )
            test_t2g_path = os.path.join(
                tempfile.mkdtemp(), 'transcripts_to_genes.txt'
            )
            with open(test_index_path, 'w') as index, open(test_t2g_path,
                                                           'w') as t2g:
                index.write('INDEX')
                t2g.write('T2G')
            test_tar_path = os.path.join(
                tempfile.gettempdir(), '{}.tar.gz'.format(uuid.uuid4())
            )
            with tarfile.open(test_tar_path, 'w:gz') as f:
                f.add(
                    test_index_path, arcname=os.path.basename(test_index_path)
                )
                f.add(test_t2g_path, arcname=os.path.basename(test_t2g_path))
            download_file.return_value = test_tar_path
            with self.assertRaises(Exception):
                ref.download_reference(reference, files, temp_dir=temp_dir)
            download_file.assert_not_called()

    def test_decompress_file_text(self):
        with mock.patch('kb_python.ref.decompress_gzip') as decompress_gzip:
            temp_dir = tempfile.mkdtemp()
            self.assertEqual(
                'textfile.txt',
                ref.decompress_file('textfile.txt', temp_dir=temp_dir)
            )
            decompress_gzip.assert_not_called()

    def test_decompress_file_gzip(self):
        with mock.patch('kb_python.ref.decompress_gzip') as decompress_gzip:
            temp_dir = tempfile.mkdtemp()
            decompress_gzip.return_value = 'textfile.txt'
            self.assertEqual(
                'textfile.txt',
                ref.decompress_file('textfile.txt.gz', temp_dir=temp_dir)
            )
            decompress_gzip.assert_called_once_with(
                'textfile.txt.gz', os.path.join(temp_dir, 'textfile.txt')
            )

    def test_ref(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.create_t2g_from_gtf') as create_t2g_from_gtf,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
            mock.patch('kb_python.ref.check_chromosomes') as check_chromosomes,\
            mock.patch('kb_python.ref.generate_cdna_fasta') as generate_cdna_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = False
            temp_dir = tempfile.mkdtemp()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            sorted_fasta_path = mock.MagicMock()
            sorted_gtf_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            chromosomes = {'1', '2'}
            decompress_file.side_effect = [self.gtf_path, self.fasta_path]
            sort_fasta.return_value = sorted_fasta_path, chromosomes
            sort_gtf.return_value = sorted_gtf_path, chromosomes
            check_chromosomes.return_value = chromosomes
            generate_cdna_fasta.return_value = cdna_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_gtf.return_value = {'t2g': t2g_path}
            self.assertEqual({
                't2g': t2g_path,
                'index': index_path,
            },
                             ref.ref(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 index_path,
                                 t2g_path,
                                 temp_dir=temp_dir
                             ))
            self.assertEqual(2, decompress_file.call_count)
            decompress_file.assert_has_calls([
                call(self.gtf_path, temp_dir=temp_dir),
                call(self.fasta_path, temp_dir=temp_dir)
            ])
            create_t2g_from_gtf.assert_called_once_with(self.gtf_path, t2g_path)
            sort_fasta.assert_called_once_with(
                self.fasta_path, os.path.join(temp_dir, SORTED_FASTA_FILENAME)
            )
            sort_gtf.assert_called_once_with(
                self.gtf_path, os.path.join(temp_dir, SORTED_GTF_FILENAME)
            )
            check_chromosomes.assert_called_once_with(chromosomes, chromosomes)
            generate_cdna_fasta.assert_called_once_with(
                sorted_fasta_path,
                sorted_gtf_path,
                cdna_fasta_path,
                chromosomes=chromosomes
            )
            kallisto_index.assert_called_once_with(
                cdna_fasta_path, index_path, k=31
            )

    def test_ref_override_k(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.create_t2g_from_gtf') as create_t2g_from_gtf,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
            mock.patch('kb_python.ref.generate_cdna_fasta') as generate_cdna_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = False
            k = 999
            temp_dir = tempfile.mkdtemp()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            sorted_fasta_path = mock.MagicMock()
            sorted_gtf_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            decompress_file.side_effect = [self.gtf_path, self.fasta_path]
            sort_fasta.return_value = sorted_fasta_path
            sort_gtf.return_value = sorted_gtf_path
            generate_cdna_fasta.return_value = cdna_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_gtf.return_value = {'t2g': t2g_path}
            self.assertEqual({
                't2g': t2g_path,
                'index': index_path,
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
            self.assertEqual(2, decompress_file.call_count)
            decompress_file.assert_has_calls([
                call(self.gtf_path, temp_dir=temp_dir),
                call(self.fasta_path, temp_dir=temp_dir)
            ])
            create_t2g_from_gtf.assert_called_once_with(self.gtf_path, t2g_path)
            sort_fasta.assert_called_once_with(
                self.fasta_path, os.path.join(temp_dir, SORTED_FASTA_FILENAME)
            )
            sort_gtf.assert_called_once_with(
                self.gtf_path, os.path.join(temp_dir, SORTED_GTF_FILENAME)
            )
            generate_cdna_fasta.assert_called_once_with(
                sorted_fasta_path, sorted_gtf_path, cdna_fasta_path
            )
            kallisto_index.assert_called_once_with(
                cdna_fasta_path, index_path, k=k
            )

    def test_ref_exists(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.create_t2g_from_gtf') as create_t2g_from_gtf,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
            mock.patch('kb_python.ref.check_chromosomes') as check_chromosomes,\
            mock.patch('kb_python.ref.generate_cdna_fasta') as generate_cdna_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = True
            cdna_fasta_path = mock.MagicMock()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            temp_dir = mock.MagicMock()
            decompress_file.return_value = self.gtf_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_gtf.return_value = {'t2g': t2g_path}
            self.assertEqual({'t2g': t2g_path},
                             ref.ref(
                                 self.fasta_path,
                                 self.gtf_path,
                                 cdna_fasta_path,
                                 index_path,
                                 t2g_path,
                                 temp_dir=temp_dir
                             ))
            decompress_file.assert_called_once_with(
                self.gtf_path, temp_dir=temp_dir
            )
            create_t2g_from_gtf.assert_called_once_with(self.gtf_path, t2g_path)
            sort_fasta.assert_not_called()
            sort_gtf.assert_not_called()
            check_chromosomes.assert_not_called()
            generate_cdna_fasta.assert_not_called()
            kallisto_index.assert_not_called()

    def test_ref_overwrite(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.create_t2g_from_gtf') as create_t2g_from_gtf,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
            mock.patch('kb_python.ref.check_chromosomes') as check_chromosomes,\
            mock.patch('kb_python.ref.generate_cdna_fasta') as generate_cdna_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = True
            temp_dir = tempfile.mkdtemp()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            sorted_fasta_path = mock.MagicMock()
            sorted_gtf_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            chromosomes = {'1', '2'}
            decompress_file.side_effect = [self.gtf_path, self.fasta_path]
            sort_fasta.return_value = sorted_fasta_path, chromosomes
            sort_gtf.return_value = sorted_gtf_path, chromosomes
            check_chromosomes.return_value = chromosomes
            generate_cdna_fasta.return_value = cdna_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_gtf.return_value = {'t2g': t2g_path}
            self.assertEqual({
                't2g': t2g_path,
                'index': index_path,
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
            self.assertEqual(2, decompress_file.call_count)
            decompress_file.assert_has_calls([
                call(self.gtf_path, temp_dir=temp_dir),
                call(self.fasta_path, temp_dir=temp_dir)
            ])
            create_t2g_from_gtf.assert_called_once_with(self.gtf_path, t2g_path)
            sort_fasta.assert_called_once_with(
                self.fasta_path, os.path.join(temp_dir, SORTED_FASTA_FILENAME)
            )
            sort_gtf.assert_called_once_with(
                self.gtf_path, os.path.join(temp_dir, SORTED_GTF_FILENAME)
            )
            check_chromosomes.assert_called_once_with(chromosomes, chromosomes)
            generate_cdna_fasta.assert_called_once_with(
                sorted_fasta_path,
                sorted_gtf_path,
                cdna_fasta_path,
                chromosomes=chromosomes
            )
            kallisto_index.assert_called_once_with(
                cdna_fasta_path, index_path, k=31
            )

    def test_ref_kite_odd(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.generate_kite_fasta') as generate_kite_fasta,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = False
            temp_dir = tempfile.mkdtemp()
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
            kallisto_index.assert_called_once_with(fasta_path, index_path, k=1)

    def test_ref_kite_even(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.generate_kite_fasta') as generate_kite_fasta,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = False
            temp_dir = tempfile.mkdtemp()
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
            kallisto_index.assert_called_once_with(fasta_path, index_path, k=1)

    def test_ref_kite_override_k(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.generate_kite_fasta') as generate_kite_fasta,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = False
            k = 999
            temp_dir = tempfile.mkdtemp()
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
            kallisto_index.assert_called_once_with(fasta_path, index_path, k=k)

    def test_ref_kite_doesnt_overwrite(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.generate_kite_fasta') as generate_kite_fasta,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = True
            temp_dir = tempfile.mkdtemp()
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
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = True
            temp_dir = tempfile.mkdtemp()
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
            kallisto_index.assert_called_once_with(fasta_path, index_path, k=1)

    def test_ref_lamanno(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
            mock.patch('kb_python.ref.check_chromosomes') as check_chromosomes,\
            mock.patch('kb_python.ref.generate_cdna_fasta') as generate_cdna_fasta,\
            mock.patch('kb_python.ref.generate_intron_fasta') as generate_intron_fasta,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = False
            temp_dir = tempfile.mkdtemp()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            sorted_fasta_path = mock.MagicMock()
            sorted_gtf_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            intron_fasta_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            chromosomes = {'1', '2'}
            decompress_file.side_effect = [self.fasta_path, self.gtf_path]
            sort_fasta.return_value = sorted_fasta_path, chromosomes
            sort_gtf.return_value = sorted_gtf_path, chromosomes
            check_chromosomes.return_value = chromosomes
            generate_cdna_fasta.return_value = cdna_fasta_path
            generate_intron_fasta.return_value = intron_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': cdna_t2c_path
            }, {
                't2c': intron_t2c_path
            }]
            concatenate_files.return_value = combined_path
            self.assertEqual({
                't2g': t2g_path,
                'cdna_fasta': cdna_fasta_path,
                'intron_fasta': intron_fasta_path,
                'cdna_t2c': cdna_t2c_path,
                'intron_t2c': intron_t2c_path,
                'index': index_path,
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
                                 temp_dir=temp_dir
                             ))
            self.assertEqual(2, decompress_file.call_count)
            decompress_file.assert_has_calls([
                call(self.fasta_path, temp_dir=temp_dir),
                call(self.gtf_path, temp_dir=temp_dir),
            ])
            create_t2g_from_fasta.assert_called_once_with(
                combined_path, t2g_path
            )
            sort_fasta.assert_called_once_with(
                self.fasta_path, os.path.join(temp_dir, SORTED_FASTA_FILENAME)
            )
            sort_gtf.assert_called_once_with(
                self.gtf_path, os.path.join(temp_dir, SORTED_GTF_FILENAME)
            )
            check_chromosomes.assert_called_once_with(chromosomes, chromosomes)
            generate_cdna_fasta.assert_called_once_with(
                sorted_fasta_path,
                sorted_gtf_path,
                cdna_fasta_path,
                chromosomes=chromosomes
            )
            generate_intron_fasta.assert_called_once_with(
                sorted_fasta_path,
                sorted_gtf_path,
                intron_fasta_path,
                chromosomes=chromosomes
            )
            self.assertEqual(2, create_t2c.call_count)
            create_t2c.assert_has_calls([
                call(cdna_fasta_path, cdna_t2c_path),
                call(intron_fasta_path, intron_t2c_path)
            ])
            kallisto_index.assert_called_once_with(
                combined_path, index_path, k=31
            )

    def test_ref_lamanno_override_k(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
            mock.patch('kb_python.ref.generate_cdna_fasta') as generate_cdna_fasta,\
            mock.patch('kb_python.ref.generate_intron_fasta') as generate_intron_fasta,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = False
            k = 999
            temp_dir = tempfile.mkdtemp()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            sorted_fasta_path = mock.MagicMock()
            sorted_gtf_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            intron_fasta_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            decompress_file.side_effect = [self.fasta_path, self.gtf_path]
            sort_fasta.return_value = sorted_fasta_path
            sort_gtf.return_value = sorted_gtf_path
            generate_cdna_fasta.return_value = cdna_fasta_path
            generate_intron_fasta.return_value = intron_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': cdna_t2c_path
            }, {
                't2c': intron_t2c_path
            }]
            concatenate_files.return_value = combined_path
            self.assertEqual({
                't2g': t2g_path,
                'cdna_fasta': cdna_fasta_path,
                'intron_fasta': intron_fasta_path,
                'cdna_t2c': cdna_t2c_path,
                'intron_t2c': intron_t2c_path,
                'index': index_path,
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
                                 k=k,
                                 temp_dir=temp_dir
                             ))
            self.assertEqual(2, decompress_file.call_count)
            decompress_file.assert_has_calls([
                call(self.fasta_path, temp_dir=temp_dir),
                call(self.gtf_path, temp_dir=temp_dir),
            ])
            create_t2g_from_fasta.assert_called_once_with(
                combined_path, t2g_path
            )
            sort_fasta.assert_called_once_with(
                self.fasta_path, os.path.join(temp_dir, SORTED_FASTA_FILENAME)
            )
            sort_gtf.assert_called_once_with(
                self.gtf_path, os.path.join(temp_dir, SORTED_GTF_FILENAME)
            )
            generate_cdna_fasta.assert_called_once_with(
                sorted_fasta_path, sorted_gtf_path, cdna_fasta_path
            )
            generate_intron_fasta.assert_called_once_with(
                sorted_fasta_path, sorted_gtf_path, intron_fasta_path
            )
            self.assertEqual(2, create_t2c.call_count)
            create_t2c.assert_has_calls([
                call(cdna_fasta_path, cdna_t2c_path),
                call(intron_fasta_path, intron_t2c_path)
            ])
            kallisto_index.assert_called_once_with(
                combined_path, index_path, k=k
            )

    def test_ref_lamanno_exists(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
            mock.patch('kb_python.ref.check_chromosomes') as check_chromosomes,\
            mock.patch('kb_python.ref.generate_cdna_fasta') as generate_cdna_fasta,\
            mock.patch('kb_python.ref.generate_intron_fasta') as generate_intron_fasta,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = True
            temp_dir = tempfile.mkdtemp()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            intron_fasta_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            self.assertEqual({},
                             ref.ref_lamanno(
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
            decompress_file.assert_not_called()
            create_t2g_from_fasta.assert_not_called()
            sort_fasta.assert_not_called()
            sort_gtf.assert_not_called()
            check_chromosomes.assert_not_called()
            generate_cdna_fasta.assert_not_called()
            generate_intron_fasta.assert_not_called()
            create_t2c.assert_not_called()
            concatenate_files.assert_not_called()
            kallisto_index.assert_not_called()

    def test_ref_lamanno_overwrite(self):
        with mock.patch('kb_python.ref.decompress_file') as decompress_file,\
            mock.patch('kb_python.ref.create_t2g_from_fasta') as create_t2g_from_fasta,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
            mock.patch('kb_python.ref.check_chromosomes') as check_chromosomes,\
            mock.patch('kb_python.ref.generate_cdna_fasta') as generate_cdna_fasta,\
            mock.patch('kb_python.ref.generate_intron_fasta') as generate_intron_fasta,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = False
            temp_dir = tempfile.mkdtemp()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            sorted_fasta_path = mock.MagicMock()
            sorted_gtf_path = mock.MagicMock()
            cdna_fasta_path = mock.MagicMock()
            intron_fasta_path = mock.MagicMock()
            cdna_t2c_path = mock.MagicMock()
            intron_t2c_path = mock.MagicMock()
            combined_path = mock.MagicMock()
            chromosomes = {'1', '2'}
            decompress_file.side_effect = [self.fasta_path, self.gtf_path]
            sort_fasta.return_value = sorted_fasta_path, chromosomes
            sort_gtf.return_value = sorted_gtf_path, chromosomes
            check_chromosomes.return_value = chromosomes
            generate_cdna_fasta.return_value = cdna_fasta_path
            generate_intron_fasta.return_value = intron_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_fasta.return_value = {'t2g': t2g_path}
            create_t2c.side_effect = [{
                't2c': cdna_t2c_path
            }, {
                't2c': intron_t2c_path
            }]
            concatenate_files.return_value = combined_path
            self.assertEqual({
                't2g': t2g_path,
                'cdna_fasta': cdna_fasta_path,
                'intron_fasta': intron_fasta_path,
                'cdna_t2c': cdna_t2c_path,
                'intron_t2c': intron_t2c_path,
                'index': index_path,
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
                                 overwrite=True
                             ))
            self.assertEqual(2, decompress_file.call_count)
            decompress_file.assert_has_calls([
                call(self.fasta_path, temp_dir=temp_dir),
                call(self.gtf_path, temp_dir=temp_dir),
            ])
            create_t2g_from_fasta.assert_called_once_with(
                combined_path, t2g_path
            )
            sort_fasta.assert_called_once_with(
                self.fasta_path, os.path.join(temp_dir, SORTED_FASTA_FILENAME)
            )
            sort_gtf.assert_called_once_with(
                self.gtf_path, os.path.join(temp_dir, SORTED_GTF_FILENAME)
            )
            check_chromosomes.assert_called_once_with(chromosomes, chromosomes)
            generate_cdna_fasta.assert_called_once_with(
                sorted_fasta_path,
                sorted_gtf_path,
                cdna_fasta_path,
                chromosomes=chromosomes
            )
            generate_intron_fasta.assert_called_once_with(
                sorted_fasta_path,
                sorted_gtf_path,
                intron_fasta_path,
                chromosomes=chromosomes
            )
            self.assertEqual(2, create_t2c.call_count)
            create_t2c.assert_has_calls([
                call(cdna_fasta_path, cdna_t2c_path),
                call(intron_fasta_path, intron_t2c_path)
            ])
            kallisto_index.assert_called_once_with(
                combined_path, index_path, k=31
            )

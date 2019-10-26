import os
import tarfile
import tempfile
import uuid
from unittest import mock, TestCase
from unittest.mock import call

import kb_python.ref as ref
from kb_python.config import INDEX_FILENAME, REFERENCES_MAPPING, T2G_FILENAME
from kb_python.constants import (
    SORTED_FASTA_FILENAME,
    SORTED_GTF_FILENAME,
)
from tests.mixins import TestMixin


class TestRef(TestMixin, TestCase):

    def test_sort_gtf(self):
        with mock.patch('kb_python.ref.GTF') as GTF:
            out_path = 'test'
            ref.sort_gtf(self.unsorted_gtf_path, out_path)
            GTF().sort.assert_called_once_with(out_path)

    def test_sort_fasta(self):
        with mock.patch('kb_python.ref.FASTA') as FASTA:
            out_path = 'test'
            ref.sort_fasta(self.unsorted_fasta_path, out_path)
            FASTA().sort.assert_called_once_with(out_path)

    def test_kallisto_index(self):
        index_path = os.path.join(
            tempfile.gettempdir(), '{}.idx'.format(uuid.uuid4())
        )
        self.assertFalse(os.path.exists(index_path))
        result = ref.kallisto_index(self.fasta_path, index_path)
        for key, path in result.items():
            self.assertTrue(os.path.exists(path))

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
        with mock.patch('kb_python.ref.urlretrieve') as urlretrieve:
            choice = 'human'
            url = REFERENCES_MAPPING[choice].url
            temp_dir = tempfile.mkdtemp()
            index_path = os.path.join(
                tempfile.mkdtemp(), '{}.idx'.format(uuid.uuid4())
            )
            t2g_path = os.path.join(
                tempfile.mkdtemp(), '{}.txt'.format(uuid.uuid4())
            )
            test_index_path = os.path.join(tempfile.mkdtemp(), INDEX_FILENAME)
            test_t2g_path = os.path.join(tempfile.mkdtemp(), T2G_FILENAME)
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
            urlretrieve.return_value = test_tar_path, None
            self.assertEqual({
                'index': index_path,
                't2g': t2g_path
            },
                             ref.download_reference(
                                 choice,
                                 index_path,
                                 t2g_path,
                                 temp_dir=temp_dir
                             ))
            urlretrieve.assert_called_once_with(
                url, os.path.join(temp_dir, os.path.basename(url))
            )
            with open(index_path, 'r') as index, open(t2g_path, 'r') as t2g:
                self.assertEqual('INDEX', index.read())
                self.assertEqual('T2G', t2g.read())

    def test_ref(self):
        with mock.patch('kb_python.ref.create_t2g_from_gtf') as create_t2g_from_gtf,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
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
                                 temp_dir=temp_dir
                             ))
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
            kallisto_index.assert_called_once_with(cdna_fasta_path, index_path)

    def test_ref_exists(self):
        with mock.patch('kb_python.ref.create_t2g_from_gtf') as create_t2g_from_gtf,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
            mock.patch('kb_python.ref.generate_cdna_fasta') as generate_cdna_fasta,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = True
            cdna_fasta_path = mock.MagicMock()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_gtf.return_value = {'t2g': t2g_path}
            self.assertEqual({'t2g': t2g_path},
                             ref.ref(
                                 self.fasta_path, self.gtf_path,
                                 cdna_fasta_path, index_path, t2g_path
                             ))
            create_t2g_from_gtf.assert_called_once_with(self.gtf_path, t2g_path)
            sort_fasta.assert_not_called()
            sort_gtf.assert_not_called()
            generate_cdna_fasta.assert_not_called()
            kallisto_index.assert_not_called()

    def test_ref_overwrite(self):
        with mock.patch('kb_python.ref.create_t2g_from_gtf') as create_t2g_from_gtf,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
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
                                 temp_dir=temp_dir,
                                 overwrite=True
                             ))
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
            kallisto_index.assert_called_once_with(cdna_fasta_path, index_path)

    def test_ref_lamanno(self):
        with mock.patch('kb_python.ref.create_t2g_from_gtf') as create_t2g_from_gtf,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
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
            sort_fasta.return_value = sorted_fasta_path
            sort_gtf.return_value = sorted_gtf_path
            generate_cdna_fasta.return_value = cdna_fasta_path
            generate_intron_fasta.return_value = intron_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_gtf.return_value = {'t2g': t2g_path}
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
            create_t2g_from_gtf.assert_called_once_with(
                self.gtf_path, t2g_path, intron=True
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
            kallisto_index.assert_called_once_with(combined_path, index_path)

    def test_ref_lamanno_exists(self):
        with mock.patch('kb_python.ref.create_t2g_from_gtf') as create_t2g_from_gtf,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
            mock.patch('kb_python.ref.generate_cdna_fasta') as generate_cdna_fasta,\
            mock.patch('kb_python.ref.generate_intron_fasta') as generate_intron_fasta,\
            mock.patch('kb_python.ref.create_t2c') as create_t2c,\
            mock.patch('kb_python.ref.concatenate_files') as concatenate_files,\
            mock.patch('kb_python.ref.kallisto_index') as kallisto_index,\
            mock.patch('kb_python.ref.os.path.exists') as exists:
            exists.return_value = True
            cdna_fasta_path = mock.MagicMock()
            index_path = mock.MagicMock()
            t2g_path = mock.MagicMock()
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_gtf.return_value = {'t2g': t2g_path}
            self.assertEqual({'t2g': t2g_path},
                             ref.ref(
                                 self.fasta_path, self.gtf_path,
                                 cdna_fasta_path, index_path, t2g_path
                             ))
            create_t2g_from_gtf.assert_called_once_with(self.gtf_path, t2g_path)
            sort_fasta.assert_not_called()
            sort_gtf.assert_not_called()
            generate_cdna_fasta.assert_not_called()
            generate_intron_fasta.assert_not_called()
            create_t2c.assert_not_called()
            concatenate_files.assert_not_called()
            kallisto_index.assert_not_called()

    def test_ref_lamanno_overwrite(self):
        with mock.patch('kb_python.ref.create_t2g_from_gtf') as create_t2g_from_gtf,\
            mock.patch('kb_python.ref.sort_fasta') as sort_fasta,\
            mock.patch('kb_python.ref.sort_gtf') as sort_gtf,\
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
            sort_fasta.return_value = sorted_fasta_path
            sort_gtf.return_value = sorted_gtf_path
            generate_cdna_fasta.return_value = cdna_fasta_path
            generate_intron_fasta.return_value = intron_fasta_path
            kallisto_index.return_value = {'index': index_path}
            create_t2g_from_gtf.return_value = {'t2g': t2g_path}
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
            create_t2g_from_gtf.assert_called_once_with(
                self.gtf_path, t2g_path, intron=True
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
            kallisto_index.assert_called_once_with(combined_path, index_path)
